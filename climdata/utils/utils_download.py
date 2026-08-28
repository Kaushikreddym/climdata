"""Google Drive helpers backing the MSWX mirror.

MSWX is distributed from a public Drive folder rather than an HTTP server, so
listing and downloading go through the Drive API. The Google client libraries are
imported optionally: :class:`~climdata.datasets.MSWX.MSWXmirror` reports their
absence itself, and importing this module must not fail for users of the other
providers.

Both functions take an already-authenticated Drive ``service`` object, built by
the caller from a service-account key.
"""

import os
import io
import warnings

try:
    from google.oauth2 import service_account  # noqa: F401 — re-exported for callers
    from googleapiclient.discovery import build  # noqa: F401 — re-exported for callers
    from googleapiclient.http import MediaIoBaseDownload
except ImportError:
    pass

warnings.filterwarnings("ignore", category=Warning)


def list_drive_files(folder_id, service):
    """List every file in a Google Drive folder, following pagination.

    The Drive API caps a listing at a few hundred entries per response, and an
    MSWX variable folder holds one file per day — tens of thousands. This walks
    ``nextPageToken`` to the end, so the result is the complete folder rather
    than its first page.

    Args:
        folder_id (str): Drive folder ID.
        service: An authenticated Drive v3 service object.

    Returns:
        list[dict]: One entry per file, each with ``id`` and ``name``. Trashed
        files are excluded.

    Raises:
        googleapiclient.errors.HttpError: If the folder does not exist or the
            service account cannot read it.
    """
    files = []
    page_token = None

    while True:
        results = service.files().list(
            q=f"'{folder_id}' in parents and trashed = false",
            fields="files(id, name), nextPageToken",
            pageToken=page_token
        ).execute()

        files.extend(results.get("files", []))
        page_token = results.get("nextPageToken", None)

        if not page_token:
            break

    return files


def download_drive_file(file_id, local_path, service):
    """Download one Drive file to a local path, streaming it in chunks.

    Parent directories are created as needed. The download is chunked rather
    than buffered whole, and progress is printed as it goes.

    A partially written file is left behind if the transfer fails, and callers
    that test for existence will then treat it as complete — delete a suspect
    file rather than re-running over it.

    Args:
        file_id (str): Drive file ID, from :func:`list_drive_files`.
        local_path (str): Destination path.
        service: An authenticated Drive v3 service object.

    Returns:
        None

    Raises:
        googleapiclient.errors.HttpError: If the file is missing or unreadable.
    """
    request = service.files().get_media(fileId=file_id)
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    with io.FileIO(local_path, 'wb') as fh:
        downloader = MediaIoBaseDownload(fh, request)

        done = False
        while not done:
            status, done = downloader.next_chunk()
            print(f"   → Download {int(status.progress() * 100)}% complete")
