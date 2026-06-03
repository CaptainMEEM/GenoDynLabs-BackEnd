"""
worker.py  --  RQ worker entry point. Run as its own Railway service:

    python worker.py

It loads the reference bundle ONCE at startup (a few MB, parsed once) and then
pulls jobs off the queue, running run_report_job for each. Each job is fully
in-memory and shares only the read-only bundle, so to handle many reports at
once you simply scale the number of worker replicas in Railway — there is no
shared mutable state and no disk contention between them.

Tuning knobs (env):
  REDIS_URL              required; the same Redis the web service enqueues to.
  RQ_QUEUES              optional; comma-separated queue names (default "reports").
"""
import os
import logging

import redis
from rq import Worker, Queue

from services.annotate import load_reference

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s %(message)s")
log = logging.getLogger("worker")

# Warm the cache before taking any job so the first report isn't slow.
_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                     "reference.pkl")
load_reference(_DATA)
log.info("reference bundle warm; ready for jobs")

if __name__ == "__main__":
    conn = redis.from_url(os.environ["REDIS_URL"])
    queues = [q.strip() for q in os.environ.get("RQ_QUEUES", "reports").split(",")]
    Worker([Queue(name, connection=conn) for name in queues],
           connection=conn).work(with_scheduler=True)
