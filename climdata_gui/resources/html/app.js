/* ======================================================================
   climdata dashboard — app.js
   Python calls:  initOverlayUI(config)  updateOverlay(payload)  clearOverlays()
   Bridge calls:  bridge.on_point_selected  on_box_selected  on_render_request
   ====================================================================== */
'use strict';

/* ── Shared state ─────────────────────────────────────────────────── */
const state = {
  map:           null,
  currentBase:   null,
  rasterOverlay: null,
  bridge:        null,
  selectionMode: 'point',   // 'point' | 'box'
  marker:        null,
  boxRect:       null,
  drawing:       false,
  startLL:       null,
  sidebarOpen:   true,
};

/* ── Basemap catalogue ────────────────────────────────────────────── */
const BASEMAPS = {
  'osm': {
    url:     'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
    options: { maxZoom: 19, attribution: '© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>' },
  },
  'carto-dark': {
    url:     'https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
    options: { maxZoom: 19, attribution: '© <a href="https://carto.com/">CartoDB</a>' },
  },
  'carto-light': {
    url:     'https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
    options: { maxZoom: 19, attribution: '© <a href="https://carto.com/">CartoDB</a>' },
  },
  'esri': {
    url:     'https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}',
    options: { maxZoom: 19, attribution: '© <a href="https://www.esri.com/">Esri</a>' },
  },
};

/* ================================================================== */
/*  Map initialisation                                                 */
/* ================================================================== */
function initMap() {
  state.map = L.map('map', {
    touchZoom:          true,
    scrollWheelZoom:    true,
    dragging:           true,
    tap:                false,
    bounceAtZoomLimits: false,
    wheelDebounceTime:  40,
    wheelPxPerZoomLevel: 60,
    zoomControl:        false,   // we add it manually at bottom-right
  }).setView([51.0, 10.0], 4);

  L.control.zoom({ position: 'bottomright' }).addTo(state.map);

  const bm = BASEMAPS['osm'];
  state.currentBase = L.tileLayer(bm.url, bm.options).addTo(state.map);

  L.Marker.prototype.options.icon = L.icon({
    iconUrl:       'marker-icon.png',
    iconRetinaUrl: 'marker-icon-2x.png',
    shadowUrl:     'marker-shadow.png',
    iconSize:  [25, 41], iconAnchor:   [12, 41],
    popupAnchor: [1, -34], shadowSize: [41, 41],
  });
}

/* ================================================================== */
/*  Sidebar & card collapse                                            */
/* ================================================================== */
function initSidebar() {
  const sidebar = document.getElementById('sidebar');
  const toggle  = document.getElementById('sidebarToggle');

  toggle.addEventListener('click', () => {
    state.sidebarOpen = !state.sidebarOpen;
    sidebar.classList.toggle('collapsed', !state.sidebarOpen);
    toggle.classList.toggle('collapsed', !state.sidebarOpen);
    toggle.textContent = state.sidebarOpen ? '◀' : '▶';
    document.body.classList.toggle('sidebar-collapsed', !state.sidebarOpen);
  });

  /* Card collapse — click anywhere on the header row */
  document.querySelectorAll('.card-header').forEach(header => {
    const chevron = header.querySelector('.card-chevron');
    const bodyId  = 'body-' + header.dataset.card;

    header.addEventListener('click', () => {
      const body = document.getElementById(bodyId);
      const closing = !body.classList.contains('closed');

      if (closing) {
        /* Lock current height then animate to 0 */
        body.style.maxHeight = body.scrollHeight + 'px';
        requestAnimationFrame(() => {
          body.classList.add('closed');
          body.style.maxHeight = '0';
        });
      } else {
        /* Open: remove closed, animate to scrollHeight, then free */
        body.classList.remove('closed');
        const target = body.scrollHeight + 'px';
        body.style.maxHeight = target;
        body.addEventListener('transitionend', () => {
          if (!body.classList.contains('closed')) body.style.maxHeight = 'none';
        }, { once: true });
      }
      chevron.classList.toggle('closed', closing);
    });
  });
}

/* ================================================================== */
/*  Tab navigation                                                     */
/* ================================================================== */
function initTabs() {
  document.body.classList.add('tab-home');

  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const tab = btn.dataset.tab;
      // Update active button
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      // Swap body tab class
      document.body.className = document.body.className
        .replace(/\btab-\S+/g, '').replace(/\s+/g, ' ').trim();
      document.body.classList.add('tab-' + tab);
      // Re-render map when switching back to Home (it may have been hidden)
      if (tab === 'download' && state.map) {
        setTimeout(() => state.map.invalidateSize(), 50);
      }
    });
  });
}

/* ================================================================== */
/*  Theme toggle                                                       */
/* ================================================================== */
function initTheme() {
  document.getElementById('themeToggle').addEventListener('click', () => {
    const html = document.documentElement;
    html.dataset.theme = html.dataset.theme === 'dark' ? 'light' : 'dark';
  });
}

/* ================================================================== */
/*  Basemap switcher                                                   */
/* ================================================================== */
function initBasemap() {
  document.querySelectorAll('.basemap-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const key = btn.dataset.basemap;
      if (!(key in BASEMAPS)) return;
      document.querySelectorAll('.basemap-btn').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      if (state.currentBase) state.map.removeLayer(state.currentBase);
      const bm = BASEMAPS[key];
      state.currentBase = L.tileLayer(bm.url, bm.options).addTo(state.map);
      if (state.rasterOverlay) state.rasterOverlay.bringToFront();
    });
  });
}

/* ================================================================== */
/*  AOI mode / clear controls                                          */
/* ================================================================== */
function initControls() {
  const modeBtn = document.getElementById('modeBtn');
  modeBtn.addEventListener('click', () => {
    const entering = state.selectionMode !== 'box';
    state.selectionMode = entering ? 'box' : 'point';
    modeBtn.dataset.active    = entering.toString();
    modeBtn.textContent       = entering ? '📍 Point' : '⬜ Box';
    state.map.getContainer().style.cursor = entering ? 'crosshair' : '';
  });

  document.getElementById('clearBtn').addEventListener('click', () => {
    clearOverlays();
    if (state.boxRect) { state.map.removeLayer(state.boxRect); state.boxRect = null; }
    if (state.marker)  { state.map.removeLayer(state.marker);  state.marker  = null; }
    setCoordText('Click map to select a point');
  });
}

/* ================================================================== */
/*  Overlay panel controls (inside sidebar)                           */
/* ================================================================== */
function initOverlayControls() {
  document.getElementById('oc-var').addEventListener('change', () => {
    document.getElementById('oc-vmin').value = '';
    document.getElementById('oc-vmax').value = '';
    requestRender();
  });

  document.getElementById('oc-cmap').addEventListener('change', requestRender);
  document.getElementById('oc-vmin').addEventListener('change', requestRender);
  document.getElementById('oc-vmax').addEventListener('change', requestRender);

  document.getElementById('oc-opacity').addEventListener('input', function () {
    const v = parseFloat(this.value);
    document.getElementById('opacity-badge').textContent = v.toFixed(2);
    if (state.rasterOverlay) state.rasterOverlay.setOpacity(v);
  });
}

/* ================================================================== */
/*  Map interaction events                                             */
/* ================================================================== */
function initMapEvents() {
  const map = state.map;

  /* ── Point selection ── */
  map.on('click', e => {
    if (state.selectionMode !== 'point') return;
    const { lat, lng } = e.latlng;
    if (state.marker) map.removeLayer(state.marker);
    state.marker = L.marker([lat, lng]).addTo(map);
    setCoordText(`Point  lat=${lat.toFixed(4)}  lon=${lng.toFixed(4)}`);
    updateAoiStatus(`Point  ${lat.toFixed(4)}, ${lng.toFixed(4)}`);
    if (state.bridge) state.bridge.on_point_selected(lat, lng);
  });

  /* ── Box rubber-band ── */
  map.on('mousedown', e => {
    if (state.selectionMode !== 'box' || e.originalEvent.button !== 0) return;
    state.drawing = true;
    state.startLL = e.latlng;
    if (state.boxRect) map.removeLayer(state.boxRect);
    state.boxRect = L.rectangle([state.startLL, state.startLL], {
      color: '#4f91ff', weight: 2,
      fillColor: '#4f91ff', fillOpacity: 0.12,
      dashArray: '6 3',
    }).addTo(map);
    map.dragging.disable();
    L.DomEvent.preventDefault(e.originalEvent);
  });

  map.on('mousemove', e => {
    if (state.drawing) state.boxRect.setBounds([state.startLL, e.latlng]);
  });

  map.on('mouseup', e => {
    if (!state.drawing) return;
    state.drawing = false;
    map.dragging.enable();
    const b = state.boxRect.getBounds();
    const sw = b.getSouthWest(), ne = b.getNorthEast();
    setCoordText(
      `Box  (${sw.lat.toFixed(3)}, ${sw.lng.toFixed(3)}) → (${ne.lat.toFixed(3)}, ${ne.lng.toFixed(3)})`
    );
    updateAoiStatus(`Box  (${sw.lat.toFixed(2)}, ${sw.lng.toFixed(2)}) → (${ne.lat.toFixed(2)}, ${ne.lng.toFixed(2)})`);
    if (state.bridge) state.bridge.on_box_selected(sw.lat, ne.lat, sw.lng, ne.lng);
  });
}

/* ================================================================== */
/*  Overlay helpers  ← called by Python via page.runJavaScript()      */
/* ================================================================== */

/** Remove raster overlay and reset the overlay panel. */
function clearOverlays() {
  if (state.rasterOverlay) {
    state.map.removeLayer(state.rasterOverlay);
    state.rasterOverlay = null;
  }
  document.getElementById('overlay-controls').style.display = 'none';
  document.getElementById('no-data-msg').style.display     = 'block';
  const img = document.getElementById('oc-colorbar-img');
  img.src = ''; img.style.display = 'none';
  document.getElementById('colorbar-container').innerHTML =
    '<img id="oc-colorbar-img" src="" alt="" style="display:none">';
}

/**
 * Initialise overlay UI: populate variable list, zoom to extent.
 * config = { vars:[...], lat_min, lon_min, lat_max, lon_max }
 */
function initOverlayUI(config) {
  console.log('[climdata] initOverlayUI', JSON.stringify(config));

  const sel = document.getElementById('oc-var');
  sel.innerHTML = '';
  config.vars.forEach(v => {
    const opt = document.createElement('option');
    opt.value = v; opt.textContent = v;
    sel.appendChild(opt);
  });

  document.getElementById('no-data-msg').style.display     = 'none';
  document.getElementById('overlay-controls').style.display = 'block';

  /* Zoom to data extent; pad single-point datasets */
  const latSpan = config.lat_max - config.lat_min;
  const lonSpan = config.lon_max - config.lon_min;
  const pad     = (latSpan < 0.1 || lonSpan < 0.1) ? 1.0 : 0;
  state.map.fitBounds([
    [config.lat_min - pad, config.lon_min - pad],
    [config.lat_max + pad, config.lon_max + pad],
  ], { padding: [30, 30] });

  setRenderLoading(true);
  requestRender();
}

/**
 * Update the raster image overlay in-place (no map pan/zoom reset).
 * payload = { b64, colorbar_b64, vmin, vmax, lat_min, lon_min, lat_max, lon_max }
 */
function updateOverlay(payload) {
  console.log('[climdata] updateOverlay b64_len=' + (payload.b64 ? payload.b64.length : 'null'));
  setRenderLoading(false);

  const dataUrl = 'data:image/png;base64,' + payload.b64;
  const bounds  = [[payload.lat_min, payload.lon_min],
                   [payload.lat_max, payload.lon_max]];
  const opacity = parseFloat(document.getElementById('oc-opacity').value);

  if (state.rasterOverlay) {
    state.rasterOverlay.setUrl(dataUrl);
    state.rasterOverlay.setBounds(bounds);
  } else {
    state.rasterOverlay = L.imageOverlay(dataUrl, bounds,
      { opacity, interactive: false, zIndex: 200 }
    ).addTo(state.map);
  }

  /* Populate vmin/vmax only when the user hasn't typed a custom value */
  if (payload.vmin !== undefined) {
    const vminEl = document.getElementById('oc-vmin');
    const vmaxEl = document.getElementById('oc-vmax');
    if (vminEl.value === '') vminEl.value = parseFloat(payload.vmin.toFixed(4));
    if (vmaxEl.value === '') vmaxEl.value = parseFloat(payload.vmax.toFixed(4));
  }

  /* Colorbar */
  if (payload.colorbar_b64) {
    const container = document.getElementById('colorbar-container');
    container.innerHTML = '<img id="oc-colorbar-img" src="" alt="">';
    const img = document.getElementById('oc-colorbar-img');
    img.src   = 'data:image/png;base64,' + payload.colorbar_b64;
    img.style.cssText = 'width:100%; display:block; border-radius:6px;';
  }
}

/**
 * Request Python to render the currently selected variable + colormap.
 * Python responds by calling updateOverlay().
 */
function requestRender() {
  if (!state.bridge || !state.bridge.on_render_request) {
    console.warn('[climdata] requestRender: bridge not ready');
    return;
  }
  const varName  = document.getElementById('oc-var').value;
  const colormap = document.getElementById('oc-cmap').value;
  const vminStr  = document.getElementById('oc-vmin').value;
  const vmaxStr  = document.getElementById('oc-vmax').value;

  console.log(`[climdata] requestRender var=${varName} cmap=${colormap}` +
              ` vmin=${vminStr || 'auto'} vmax=${vmaxStr || 'auto'}`);

  if (varName) {
    setRenderLoading(true);
    state.bridge.on_render_request(
      varName, colormap,
      vminStr !== '' ? vminStr : 'auto',
      vmaxStr !== '' ? vmaxStr : 'auto',
    );
  }
}

/* ================================================================== */
/*  Helper utilities                                                   */
/* ================================================================== */

function setCoordText(text) {
  document.getElementById('coord-display').textContent = text;
}

function updateAoiStatus(text) {
  const el = document.getElementById('dl-aoi-status');
  if (!el) return;
  el.textContent = text;
  el.classList.toggle('ready', !!text && text !== 'No AOI selected');
}

function setRenderLoading(isLoading) {
  const container = document.getElementById('colorbar-container');
  if (isLoading) {
    container.innerHTML =
      '<div class="render-loading">' +
      '<span class="loading-ring"></span>Rendering…</div>';
  }
}

/* ================================================================== */
/*  Toolbar (dataset, dates, format, data-dir, run/plot/adv buttons)  */
/* ================================================================== */
function initToolbar() {
  document.getElementById('tb-dataset').addEventListener('change', e => {
    if (state.bridge) state.bridge.on_dataset_changed(e.target.value);
  });

  function sendDates() {
    const start = document.getElementById('tb-start').value;
    const end   = document.getElementById('tb-end').value;
    if (start && end && state.bridge) state.bridge.on_dates_changed(start, end);
  }
  document.getElementById('tb-start').addEventListener('change', sendDates);
  document.getElementById('tb-end').addEventListener('change', sendDates);

  document.getElementById('tb-format').addEventListener('change', e => {
    if (state.bridge) state.bridge.on_format_changed(e.target.value);
  });

  document.getElementById('tb-dir-btn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_data_dir_browse();
  });

  document.getElementById('tb-adv-btn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_advanced_clicked();
  });

  document.getElementById('tb-run-btn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_run_clicked();
  });

  document.getElementById('tb-plot-btn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_plot_clicked();
  });
}

/* ================================================================== */
/*  Log panel                                                          */
/* ================================================================== */
function initLog() {
  document.getElementById('log-header').addEventListener('click', e => {
    if (e.target.id === 'log-clear-btn' || e.target.id === 'log-toggle-btn') return;
    toggleLog();
  });
  document.getElementById('log-toggle-btn').addEventListener('click', e => {
    e.stopPropagation();
    toggleLog();
  });
  document.getElementById('log-clear-btn').addEventListener('click', e => {
    e.stopPropagation();
    document.getElementById('log-output').innerHTML = '';
    setLogIndicator('idle');
  });
}

function toggleLog(forceOpen) {
  const open = forceOpen !== undefined ? forceOpen : !document.body.classList.contains('log-open');
  document.body.classList.toggle('log-open', open);
  const btn = document.getElementById('log-toggle-btn');
  if (btn) btn.textContent = open ? '▴' : '▾';
}

function setLogIndicator(s) {
  const el = document.getElementById('log-indicator');
  el.className = 'log-indicator' + (s !== 'idle' ? ' ' + s : '');
}

function appendLog(msg) {
  const output = document.getElementById('log-output');
  const line   = document.createElement('div');
  line.className = 'log-line';
  if (/error|fail|exception/i.test(msg))              { line.classList.add('error');   setLogIndicator('error'); }
  else if (/warn/i.test(msg))                          { line.classList.add('warn'); }
  else if (/complet|success|ready|✓/i.test(msg))      { line.classList.add('success'); setLogIndicator('active'); }
  else                                                 { setLogIndicator('active'); }
  line.textContent = msg;
  output.appendChild(line);
  output.scrollTop = output.scrollHeight;
  if (!document.body.classList.contains('log-open')) toggleLog(true);
}

/* ── State sync helpers ─────────────────────────────────────────── */
function setRunState(running) {
  const btn = document.getElementById('tb-run-btn');
  if (!btn) return;
  btn.disabled    = running;
  btn.textContent = running ? '⏳ Running…' : '▶ Run';
  btn.classList.toggle('running', running);
}

function setPlotEnabled(enabled) {
  const btn = document.getElementById('tb-plot-btn');
  if (btn) btn.disabled = !enabled;
}

function setDataDir(path) {
  const el = document.getElementById('tb-data-dir-display');
  if (!el) return;
  el.textContent = path || 'default';
  el.title = path || 'Default from config';
}

function populateDatasets(datasets) {
  const selIds = ['tb-dataset', 'basd-ref-dataset', 'basd-target-dataset', 'cmp-a-dataset', 'cmp-b-dataset'];
  selIds.forEach(id => {
    const sel = document.getElementById(id);
    if (!sel) return;
    const cur = sel.value;
    sel.innerHTML = '';
    datasets.forEach(name => {
      const opt = document.createElement('option');
      opt.value = name; opt.textContent = name;
      sel.appendChild(opt);
    });
    if (cur && datasets.includes(cur)) sel.value = cur;
  });
}

function syncToolbar(payload) {
  if (payload.dataset) {
    const sel = document.getElementById('tb-dataset');
    if (sel && sel.querySelector(`option[value="${payload.dataset}"]`)) sel.value = payload.dataset;
  }
  if (payload.start) document.getElementById('tb-start').value = payload.start;
  if (payload.end)   document.getElementById('tb-end').value   = payload.end;
  if (payload.data_dir !== undefined) setDataDir(payload.data_dir);
}

function onDashboardReady(payload) {
  populateDatasets(payload.datasets || []);
  syncToolbar(payload);
  // Confirm the actual displayed value back to Python so AppState stays in sync
  // even if this call races with an in-flight on_dataset_changed message.
  const sel = document.getElementById('tb-dataset');
  if (sel && state.bridge) state.bridge.on_dataset_changed(sel.value);
}

/* ================================================================== */
/*  QWebChannel bridge                                                 */
/* ================================================================== */
function initBridge() {
  if (typeof QWebChannel !== 'undefined') {
    new QWebChannel(qt.webChannelTransport, channel => {
      state.bridge = channel.objects.bridge;
      console.log('[climdata] QWebChannel bridge ready');

      /* Connect Python → JS signals */
      state.bridge.log_message.connect(appendLog);
      state.bridge.run_state_changed.connect(setRunState);
      state.bridge.plot_enabled.connect(setPlotEnabled);
      state.bridge.data_dir_chosen.connect(setDataDir);

      /* Notify Python that the page is ready (deferred on Python side) */
      state.bridge.on_page_ready();
    });
  } else {
    console.warn('[climdata] QWebChannel unavailable — running outside Qt?');
  }
}

/* ================================================================== */
/*  Bootstrap                                                          */
/* ================================================================== */
(function bootstrap() {
  initMap();
  initTabs();
  initSidebar();
  initTheme();
  initBasemap();
  initControls();
  initOverlayControls();
  initMapEvents();
  initToolbar();
  initLog();
  initBridge();
})();

