"""RQ worker entry point. Run as a separate Railway service:  python worker.py
Loads the 2 MB reference once at startup, then serves jobs from the queue."""
import os, redis
from rq import Worker, Queue, Connection
from services.annotate import load_reference   # warm the cache once

load_reference(os.path.join(os.path.dirname(__file__), "services", "data", "reference.pkl"))

if __name__ == "__main__":
    conn = redis.from_url(os.environ["REDIS_URL"])
    with Connection(conn):
        Worker([Queue("reports")]).work()
