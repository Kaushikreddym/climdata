import os
import io
from datetime import datetime
import warnings

import pandas as pd
import numpy as np
from omegaconf import DictConfig

from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

warnings.filterwarnings("ignore", category=Warning)

def list_drive_files(folder_id, service):
    """
    List all files in a Google Drive folder, handling pagination.
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
    """
    Download a single file from Drive to a local path.
    """
    request = service.files().get_media(fileId=file_id)
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    with io.FileIO(local_path, 'wb') as fh:
        downloader = MediaIoBaseDownload(fh, request)

        done = False
        while not done:
            status, done = downloader.next_chunk()
            print(f"   → Download {int(status.progress() * 100)}% complete")
