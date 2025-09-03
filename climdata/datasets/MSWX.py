from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
import datetime
class MSWXmirror:
    def __init__(self,cfg):
        self.cfg = cfg
        self.provider = cfg.dataset.lower()
        self.parameter_key = cfg.weather.parameter
        self.lat = cfg.location.lat
        self.lon = cfg.location.lon
        self.start_date = datetime.fromisoformat(cfg.time_range.start_date)
        self.end_date = datetime.fromisoformat(cfg.time_range.end_date)
        self.output_dir = cfg.data_dir

        provider_cfg = cfg.mappings[self.provider]
        self.param_info = provider_cfg['variables'][self.parameter_key]
        self.folder_id = self.param_info['folder_id']
        self.units = self.param_info.get("units", "")
        self.service = self._build_drive_service(provider_cfg.params.google_service_account)

    def _list_drive_files(folder_id, service):
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
    def _download_drive_file(file_id, local_path, service):
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
    def fetch():
        expected_files = []
        current = self.start_date
        while current <= self.end_date:
            doy = current.timetuple().tm_yday
            basename = f"{current.year}{doy:03d}.nc"
            expected_files.append(basename)
            current += timedelta(days=1)

        output_dir = var_cfg.data_dir
        local_files = []
        missing_files = []

        for basename in expected_files:
            local_path = os.path.join(output_dir, provider, parameter_key, basename)
            if os.path.exists(local_path):
                local_files.append(basename)
            else:
                missing_files.append(basename)

        if not missing_files:
            print(f"✅ All {len(expected_files)} files already exist locally. No download needed.")
            return local_files

        print(f"📂 {len(local_files)} exist, {len(missing_files)} missing — fetching from Drive...")

        # === 2) Connect to Drive ===
        SCOPES = ['https://www.googleapis.com/auth/drive.readonly']
        creds = service_account.Credentials.from_service_account_file(
            param_mapping[provider].params.google_service_account, scopes=SCOPES
        )
        service = build('drive', 'v3', credentials=creds)

        # === 3) List all Drive files ===
        drive_files = list_drive_files(folder_id, service)
        valid_filenames = set(missing_files)

        files_to_download = [f for f in drive_files if f['name'] in valid_filenames]

        if not files_to_download:
            print(f"⚠️ None of the missing files found in Drive. Check folder & date range.")
            return local_files

        # === 4) Download missing ===
        for file in files_to_download:
            filename = file['name']
            local_path = os.path.join(output_dir, provider, parameter_key, filename)
            print(f"⬇️ Downloading {filename} ...")
            download_drive_file(file['id'], local_path, service)
            local_files.append(filename)

        return local_files

def extract_ts_MSWX(cfg: DictConfig):
    parameter = cfg.weather.parameter
    param_mapping = cfg.mappings
    provider = cfg.dataset.lower()
    parameter_key = cfg.weather.parameter
    # Validate provider and parameter

    param_info = param_mapping[provider]['variables'][parameter_key]

    base_dir = cfg.data_dir

    target_lat = cfg.location.lat
    target_lon = cfg.location.lon

    start_date = pd.to_datetime(cfg.time_range.start_date)
    end_date = pd.to_datetime(cfg.time_range.end_date)

    # === 1) Rebuild exact basenames ===
    current = start_date
    basenames = []
    while current <= end_date:
        doy = current.timetuple().tm_yday
        basename = f"{current.year}{doy:03d}.nc"
        basenames.append(basename)
        current += timedelta(days=1)

    # === 2) Process only those files ===
    ts_list = []
    missing = []

    for basename in basenames:
        file_path = os.path.join(base_dir, provider, parameter, basename)

        if not os.path.exists(file_path):
            missing.append(basename)
            continue

        print(f"📂 Opening: {file_path}")
        ds = xr.open_dataset(file_path)

        time_name = [x for x in ds.coords if "time" in x.lower()][0]
        data_var = [v for v in ds.data_vars][0]

        ts = ds[data_var].sel(
            lat=target_lat,
            lon=target_lon,
            method='nearest'
        )

        df = ts.to_dataframe().reset_index()[[time_name, data_var]]
        ts_list.append(df)

    if missing:
        print(f"⚠️ Warning: {len(missing)} files were missing and skipped:")
        for m in missing:
            print(f"   - {m}")

    if not ts_list:
        raise RuntimeError("❌ No valid files were found. Cannot extract time series.")

    # === 3) Combine and slice (for safety) ===
    ts_all = pd.concat(ts_list).sort_values(by=time_name).reset_index(drop=True)

    ts_all[time_name] = pd.to_datetime(ts_all[time_name])
    ts_all = ts_all[
        (ts_all[time_name] >= start_date) &
        (ts_all[time_name] <= end_date)
    ].reset_index(drop=True)

    out_dir = hydra.utils.to_absolute_path(cfg.output.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, cfg.output.filename)

    ts_all["variable"] = param_info['name']
    ts_all["latitude"] = target_lat
    ts_all["longitude"] = target_lon
    ts_all['source'] = provider.upper()
    ts_all['units'] = ts.attrs['units']
    ts_all.rename(columns={param_info['name']: 'value'}, inplace=True)
    ts_all = ts_all[["latitude", "longitude", "time", "source", "variable", "value",'units']]
    
    ts_all.to_csv(out_path, index=False)
    print(f"✅ Saved MSWX time series to: {out_path}")

    return ts_all

