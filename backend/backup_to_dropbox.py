#!/usr/bin/env python
"""Backup Heroku PostgreSQL database to Dropbox.

This script:
1. Captures a Heroku database backup
2. Downloads it
3. Uploads to Dropbox
4. Cleans up old backups (keeps last 30 days)

Usage:
    python backend/backup_to_dropbox.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from datetime import UTC, datetime, timedelta

import dropbox
from dropbox.exceptions import ApiError
from dropbox.files import WriteMode

# Configuration
HEROKU_APP_NAME = os.environ.get("HEROKU_APP_NAME", "fpbase")
DROPBOX_ACCESS_TOKEN = os.environ.get("DROPBOX_ACCESS_TOKEN")
DROPBOX_BACKUP_PATH = os.environ.get("DROPBOX_BACKUP_PATH", "/heroku-backups/fpbase")
BACKUP_RETENTION_DAYS = int(os.environ.get("BACKUP_RETENTION_DAYS", "60"))


def log(message: str) -> None:
    """Print timestamped log message."""
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{timestamp}] {message}")


def run_command(cmd: list[str]) -> tuple[int, str, str]:
    """Run shell command and return exit code, stdout, stderr."""
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.returncode, result.stdout, result.stderr


def capture_heroku_backup() -> str | None:
    """Capture Heroku database backup and return backup ID."""
    log("Capturing Heroku database backup...")
    cmd = ["heroku", "pg:backups:capture", "DATABASE_URL", "--app", HEROKU_APP_NAME]
    exit_code, _stdout, stderr = run_command(cmd)

    if exit_code != 0:
        log(f"ERROR: Failed to capture backup: {stderr}")
        return None

    log("Backup captured successfully")
    return "latest"


def get_backup_url() -> str | None:
    """Get download URL for the latest backup."""
    log("Getting backup download URL...")
    cmd = ["heroku", "pg:backups:url", "--app", HEROKU_APP_NAME]
    exit_code, stdout, stderr = run_command(cmd)

    if exit_code != 0:
        log(f"ERROR: Failed to get backup URL: {stderr}")
        return None

    url = stdout.strip()
    log("Got backup URL")
    return url


def download_backup(url: str, local_path: str) -> bool:
    """Download backup from URL to local path."""
    log(f"Downloading backup to {local_path}...")
    cmd = ["curl", "-o", local_path, url]
    exit_code, _stdout, stderr = run_command(cmd)

    if exit_code != 0:
        log(f"ERROR: Failed to download backup: {stderr}")
        return False

    # Check if file was downloaded
    if not os.path.exists(local_path) or os.path.getsize(local_path) == 0:
        log("ERROR: Downloaded file is empty or doesn't exist")
        return False

    file_size_mb = os.path.getsize(local_path) / (1024 * 1024)
    log(f"Download complete ({file_size_mb:.2f} MB)")
    return True


def upload_to_dropbox(local_path: str, dropbox_path: str) -> bool:
    """Upload backup file to Dropbox."""
    if not DROPBOX_ACCESS_TOKEN:
        log("ERROR: DROPBOX_ACCESS_TOKEN environment variable not set")
        return False

    log(f"Uploading to Dropbox: {dropbox_path}...")

    try:
        dbx = dropbox.Dropbox(DROPBOX_ACCESS_TOKEN)

        # Read file in chunks to handle large files
        chunk_size = 4 * 1024 * 1024  # 4 MB chunks
        file_size = os.path.getsize(local_path)

        with open(local_path, "rb") as f:
            if file_size <= chunk_size:
                # Small file - upload in one request
                dbx.files_upload(f.read(), dropbox_path, mode=WriteMode("overwrite"))
            else:
                # Large file - use upload session
                upload_session_start_result = dbx.files_upload_session_start(f.read(chunk_size))
                cursor = dropbox.files.UploadSessionCursor(
                    session_id=upload_session_start_result.session_id,
                    offset=f.tell(),
                )
                commit = dropbox.files.CommitInfo(path=dropbox_path, mode=WriteMode("overwrite"))

                while f.tell() < file_size:
                    if (file_size - f.tell()) <= chunk_size:
                        dbx.files_upload_session_finish(f.read(chunk_size), cursor, commit)
                    else:
                        dbx.files_upload_session_append_v2(f.read(chunk_size), cursor)
                        cursor.offset = f.tell()

        log(f"Upload complete: {dropbox_path}")
        return True

    except ApiError as e:
        log(f"ERROR: Dropbox API error: {e}")
        return False
    except Exception as e:
        log(f"ERROR: Unexpected error during upload: {e}")
        return False


def cleanup_old_backups() -> None:
    """Delete Dropbox backups older than BACKUP_RETENTION_DAYS."""
    if not DROPBOX_ACCESS_TOKEN:
        return

    log(f"Cleaning up backups older than {BACKUP_RETENTION_DAYS} days...")

    try:
        dbx = dropbox.Dropbox(DROPBOX_ACCESS_TOKEN)
        cutoff_date = datetime.now(UTC) - timedelta(days=BACKUP_RETENTION_DAYS)

        # List all files in backup directory
        result = dbx.files_list_folder(DROPBOX_BACKUP_PATH)
        deleted_count = 0

        for entry in result.entries:
            if isinstance(entry, dropbox.files.FileMetadata):
                # Check if file is older than retention period
                # Make cutoff_date timezone-aware to match Dropbox's server_modified
                if entry.server_modified.replace(tzinfo=None) < cutoff_date.replace(tzinfo=None):
                    dbx.files_delete_v2(entry.path_display)
                    log(f"Deleted old backup: {entry.name}")
                    deleted_count += 1

        if deleted_count > 0:
            log(f"Cleaned up {deleted_count} old backup(s)")
        else:
            log("No old backups to clean up")

    except ApiError as e:
        if e.error.is_path() and e.error.get_path().is_not_found():
            log(f"Backup directory doesn't exist yet: {DROPBOX_BACKUP_PATH}")
        else:
            log(f"ERROR: Failed to clean up old backups: {e}")
    except Exception as e:
        log(f"ERROR: Unexpected error during cleanup: {e}")


def main() -> int:
    """Main backup workflow."""
    log("=" * 60)
    log("Starting Heroku PostgreSQL -> Dropbox backup")
    log("=" * 60)

    # Generate filename with timestamp
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d-%H%M%S")
    filename = f"fpbase-{timestamp}.dump"
    local_path = f"/tmp/{filename}"
    dropbox_path = f"{DROPBOX_BACKUP_PATH}/{filename}"

    try:
        # Step 1: Capture backup
        backup_id = capture_heroku_backup()
        if not backup_id:
            return 1

        # Step 2: Get download URL
        url = get_backup_url()
        if not url:
            return 1

        # Step 3: Download backup
        if not download_backup(url, local_path):
            return 1

        # Step 4: Upload to Dropbox
        if not upload_to_dropbox(local_path, dropbox_path):
            return 1

        # Step 5: Cleanup old backups
        cleanup_old_backups()

        log("=" * 60)
        log("Backup completed successfully!")
        log("=" * 60)
        return 0

    except Exception as e:
        log(f"ERROR: Unexpected error: {e}")
        return 1

    finally:
        # Clean up local file
        if os.path.exists(local_path):
            os.remove(local_path)
            log(f"Cleaned up local file: {local_path}")


if __name__ == "__main__":
    sys.exit(main())
