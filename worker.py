"""
worker.py  --  RQ worker entry point. Run as its own Railway service:

    python worker.py

Set WORKER_PROCESSES to how many report renders this replica should run at once
(each render is ~1 CPU core for ~20s, and holds the reference bundle in RAM, so
size it to the replica's cores AND memory). On an 8 vCPU / 8 GB Railway replica,
6 is a good default (~1.2 GB for the bundles + render headroom, 2 cores spare).

Scale horizontally by raising the REPLICA COUNT of this service in Railway. All
replicas point at the same REDIS_URL and pull DISTINCT jobs from the one queue
(RQ moves each job atomically to exactly one worker), so 5 replicas x 6 = 30
reports run concurrently with no duplicate work and no coordination code.

Env:
  REDIS_URL          required; the same Redis the web service enqueues to.
  WORKER_PROCESSES   processes per replica (default 6). Set 1 for a single worker.
  RQ_QUEUES          comma-separated queue names (default "reports").
"""
import os
import logging

import redis

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s %(message)s")
log = logging.getLogger("worker")

QUEUES = [q.strip() for q in os.environ.get("RQ_QUEUES", "reports").split(",")]
N = int(os.environ.get("WORKER_PROCESSES", "6"))


def main():
    conn = redis.from_url(os.environ["REDIS_URL"])

    if N <= 1:
        # Single process: warm the reference cache up front so the first report
        # isn't slow, then serve the queue.
        from rq import Worker, Queue
        from services.annotate import load_reference
        _data = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "data", "reference.pkl")
        load_reference(_data)
        log.info("reference bundle warm; serving %s (1 process)", QUEUES)
        Worker([Queue(q, connection=conn) for q in QUEUES],
               connection=conn).work(with_scheduler=True)
    else:
        # Multi-process pool: spawn N worker processes on this replica. Each
        # process loads the reference bundle on its first job (~0.5s) and then
        # reuses it, using the replica's cores for parallel renders.
        from rq.worker_pool import WorkerPool
        log.info("starting WorkerPool: %d processes on %s", N, QUEUES)
        WorkerPool(QUEUES, connection=conn, num_workers=N).start()


if __name__ == "__main__":
    main()
