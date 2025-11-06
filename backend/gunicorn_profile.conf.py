"""
Gunicorn configuration for memory profiling.
Mimics Heroku environment with profiling hooks.
"""
import multiprocessing
import os
import tracemalloc
from datetime import datetime

# Match Heroku settings
workers = int(os.getenv("WEB_CONCURRENCY", 2))
worker_class = "sync"  # or "gthread"
threads = 3  # only if worker_class = "gthread"
max_requests = 500
max_requests_jitter = 50
timeout = 30
keepalive = 5
preload_app = True

# Local binding
bind = "0.0.0.0:8000"
chdir = os.path.join(os.path.dirname(__file__))

# Logging
accesslog = "-"
errorlog = "-"
loglevel = "info"


# Memory profiling hooks
def on_starting(server):
    """Called just before the master process is initialized."""
    print(f"[PROFILE] Master process starting (PID: {os.getpid()})")
    tracemalloc.start()


def post_fork(server, worker):
    """Called just after a worker has been forked."""
    worker_pid = os.getpid()
    print(f"[PROFILE] Worker {worker.pid} forked (PID: {worker_pid})")

    # Start tracemalloc in each worker
    tracemalloc.start()

    # Optional: Start memray in each worker
    # import memray
    # tracker = memray.Tracker(f"memray-worker-{worker_pid}.bin")
    # tracker.__enter__()
    # worker._memray_tracker = tracker


def pre_request(worker, req):
    """Called just before a worker processes a request."""
    # Track memory before request
    import psutil
    process = psutil.Process()
    worker._request_memory_start = process.memory_info().rss / 1024 / 1024  # MB


def post_request(worker, req, environ, resp):
    """Called after a worker processes the request."""
    import psutil
    process = psutil.Process()
    mem_end = process.memory_info().rss / 1024 / 1024  # MB
    mem_start = getattr(worker, '_request_memory_start', mem_end)
    mem_delta = mem_end - mem_start

    if mem_delta > 1:  # Only log if significant increase
        print(f"[PROFILE] Worker {worker.pid}: {req.path} used {mem_delta:.2f}MB "
              f"(total: {mem_end:.2f}MB)")


def worker_exit(server, worker):
    """Called just after a worker has been exited."""
    print(f"[PROFILE] Worker {worker.pid} exiting")

    # Print top memory allocations
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    print(f"\n[PROFILE] Top 10 memory allocations in worker {worker.pid}:")
    for stat in top_stats[:10]:
        print(f"  {stat}")

    tracemalloc.stop()

    # Close memray if used
    # if hasattr(worker, '_memray_tracker'):
    #     worker._memray_tracker.__exit__(None, None, None)


def on_exit(server):
    """Called just before the master process exits."""
    print("[PROFILE] Master process exiting")
    tracemalloc.stop()
