#!/usr/bin/env python3
"""
Load testing script for local Gunicorn profiling.
Simulates realistic traffic patterns to FPbase endpoints.

Usage:
    python load_test.py --duration 60 --concurrent 10
"""
import argparse
import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests


# Common FPbase endpoints to test
ENDPOINTS = [
    "/",
    "/api/proteins/",
    "/api/proteins/?fields=name,slug",
    "/graphql/",  # Will need POST with query
    "/proteins/",
    "/about/",
    "/table/",
]

# GraphQL queries to test
GRAPHQL_QUERIES = [
    """
    query {
        proteins(first: 10) {
            edges {
                node {
                    name
                    slug
                }
            }
        }
    }
    """,
    """
    query {
        protein(name: "EGFP") {
            name
            seq
            agg
        }
    }
    """,
]


def make_request(base_url, endpoint):
    """Make a single request and return timing/size info."""
    start = time.time()
    try:
        if endpoint == "/graphql/":
            # Random GraphQL query
            query = random.choice(GRAPHQL_QUERIES)
            response = requests.post(
                f"{base_url}{endpoint}",
                json={"query": query},
                timeout=30,
            )
        else:
            response = requests.get(f"{base_url}{endpoint}", timeout=30)

        duration = time.time() - start
        return {
            "endpoint": endpoint,
            "status": response.status_code,
            "duration": duration,
            "size": len(response.content),
        }
    except Exception as e:
        duration = time.time() - start
        return {
            "endpoint": endpoint,
            "status": "ERROR",
            "duration": duration,
            "error": str(e),
        }


def run_load_test(base_url, duration_seconds, concurrent_requests):
    """Run load test for specified duration."""
    print(f"Starting load test against {base_url}")
    print(f"Duration: {duration_seconds}s, Concurrent requests: {concurrent_requests}")
    print("-" * 80)

    end_time = time.time() + duration_seconds
    total_requests = 0
    successful_requests = 0
    total_duration = 0
    errors = []

    with ThreadPoolExecutor(max_workers=concurrent_requests) as executor:
        while time.time() < end_time:
            # Submit batch of requests
            futures = []
            for _ in range(concurrent_requests):
                endpoint = random.choice(ENDPOINTS)
                future = executor.submit(make_request, base_url, endpoint)
                futures.append(future)

            # Collect results
            for future in as_completed(futures):
                result = future.result()
                total_requests += 1

                if result["status"] == 200:
                    successful_requests += 1
                    total_duration += result["duration"]
                else:
                    errors.append(result)

                # Print periodic updates
                if total_requests % 50 == 0:
                    avg_duration = total_duration / successful_requests if successful_requests else 0
                    print(
                        f"Requests: {total_requests}, "
                        f"Success: {successful_requests}, "
                        f"Avg time: {avg_duration:.3f}s"
                    )

            # Small delay between batches
            time.sleep(0.5)

    # Final summary
    print("\n" + "=" * 80)
    print("LOAD TEST SUMMARY")
    print("=" * 80)
    print(f"Total requests: {total_requests}")
    print(f"Successful: {successful_requests}")
    print(f"Failed: {len(errors)}")
    print(f"Success rate: {successful_requests / total_requests * 100:.1f}%")
    if successful_requests:
        print(f"Average response time: {total_duration / successful_requests:.3f}s")

    if errors:
        print(f"\nErrors encountered: {len(errors)}")
        for error in errors[:5]:  # Show first 5 errors
            print(f"  {error['endpoint']}: {error.get('error', error['status'])}")


def main():
    parser = argparse.ArgumentParser(description="Load test FPbase locally")
    parser.add_argument(
        "--url",
        default="http://localhost:8000",
        help="Base URL (default: http://localhost:8000)",
    )
    parser.add_argument(
        "--duration",
        type=int,
        default=60,
        help="Test duration in seconds (default: 60)",
    )
    parser.add_argument(
        "--concurrent",
        type=int,
        default=10,
        help="Concurrent requests (default: 10)",
    )

    args = parser.parse_args()

    try:
        run_load_test(args.url, args.duration, args.concurrent)
    except KeyboardInterrupt:
        print("\n\nLoad test interrupted by user")


if __name__ == "__main__":
    main()
