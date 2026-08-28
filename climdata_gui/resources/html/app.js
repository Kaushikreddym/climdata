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
  cmipDatasets:  ['cmip'],   // datasets driven by a CMIP6 model + experiment
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
/* Tabs that show the map and share the area-of-interest selection. */
const MAP_TABS = ['download', 'basd', 'comparison'];

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
        .replace(/\btab-\S+/g, '').replace(/\bmap-tab\b/g, '')
        .replace(/\s+/g, ' ').trim();
      document.body.classList.add('tab-' + tab);
      // The map-backed tabs share one map, AOI and coordinate readout
      if (MAP_TABS.includes(tab)) {
        document.body.classList.add('map-tab');
        // The map was hidden (or a different width) — recompute its size
        if (state.map) setTimeout(() => state.map.invalidateSize(), 50);
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
/*  Restart — interrupt every running job and start a fresh process    */
/* ================================================================== */
function initRestart() {
  document.getElementById('restartBtn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_restart_requested();
  });
}

/** Called by Python once the user confirms — the window is about to go. */
function setRestarting(pending) {
  const btn = document.getElementById('restartBtn');
  if (!btn) return;
  btn.disabled    = pending;
  btn.textContent = pending ? '⟳ Restarting…' : '⟳ Restart';
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
    updateAoiStatus('No AOI selected');
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

/** Mirror the current AOI into every tab that consumes it. */
function updateAoiStatus(text) {
  ['dl-aoi-status', 'basd-aoi-status', 'cmp-aoi-status'].forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;
    el.textContent = text;
    el.classList.toggle('ready', !!text && text !== 'No AOI selected');
  });
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
    updateCmipVisibility('download', e.target.value);
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
/*  BASD panel                                                         */
/* ================================================================== */
function initBasd() {
  document.getElementById('basd-out-dir-btn').addEventListener('click', () => {
    if (state.bridge) state.bridge.on_data_dir_browse();
  });

  document.getElementById('basd-run-btn').addEventListener('click', () => {
    if (!state.bridge) return;
    const cmip = cmipSelection('basd-target');
    state.bridge.on_basd_run(JSON.stringify({
      ref_dataset:    document.getElementById('basd-ref-dataset').value,
      ref_start:      document.getElementById('basd-ref-start').value,
      ref_end:        document.getElementById('basd-ref-end').value,
      target_dataset: document.getElementById('basd-target-dataset').value,
      target_start:   document.getElementById('basd-target-start').value,
      target_end:     document.getElementById('basd-target-end').value,
      method:         document.getElementById('basd-method').value,
      variable:       document.getElementById('basd-variable').value,
      out_format:     document.getElementById('basd-format').value,
      experiment_id:  cmip.experiment,
      source_id:      cmip.model,
    }));
  });
}

function setBasdRunState(running) {
  const btn = document.getElementById('basd-run-btn');
  if (!btn) return;
  btn.disabled    = running;
  btn.textContent = running ? '⏳ Running…' : '▶ Run BASD';
  btn.classList.toggle('running', running);
}

/**
 * Show a finished BASD run  ← called by Python.
 * payload = { method, variable, filename, metrics:{before,after}, notes:[] }
 */
function setBasdResult(payload) {
  const box  = document.getElementById('basd-results');
  const body = document.getElementById('basd-result-body');
  if (!box || !body) return;
  body.innerHTML = '';

  const m = payload.metrics || {};
  if (m.before && m.after) {
    body.appendChild(metricsTable(
      'Validation over the reference period',
      ['Metric', 'Raw model', 'Corrected'],
      [
        ['Bias',        m.before.model_bias, m.after.model_bias, 'abs'],
        ['RMSE',        m.before.rmse,       m.after.rmse,       'abs'],
        ['MAE',         m.before.mae,        m.after.mae,        'abs'],
        ['Correlation', m.before.correlation, m.after.correlation, 'high'],
      ]));
  }

  const summary = document.createElement('div');
  summary.innerHTML =
    `<strong>${escapeHtml(String(payload.method || '').toUpperCase())}</strong> ` +
    `applied to <strong>${escapeHtml(payload.variable || '')}</strong>.`;
  body.appendChild(summary);

  if (payload.filename) {
    const f = document.createElement('div');
    f.className = 'result-file';
    f.textContent = '💾 ' + payload.filename;
    body.appendChild(f);
  }
  (payload.notes || []).forEach(text => {
    const n = document.createElement('div');
    n.className = 'result-note';
    n.textContent = '• ' + text;
    body.appendChild(n);
  });

  box.style.display = '';
}

/* ================================================================== */
/*  Comparison panel                                                   */
/* ================================================================== */
function initComparison() {
  document.getElementById('cmp-run-btn').addEventListener('click', () => {
    if (!state.bridge) return;
    const a = cmipSelection('cmp-a');
    const b = cmipSelection('cmp-b');
    state.bridge.on_comparison_run(JSON.stringify({
      a_dataset:  document.getElementById('cmp-a-dataset').value,
      a_variable: document.getElementById('cmp-a-var').value,
      a_start:    document.getElementById('cmp-a-start').value,
      a_end:      document.getElementById('cmp-a-end').value,
      b_dataset:  document.getElementById('cmp-b-dataset').value,
      b_variable: document.getElementById('cmp-b-var').value,
      b_start:    document.getElementById('cmp-b-start').value,
      b_end:      document.getElementById('cmp-b-end').value,
      a_experiment: a.experiment, a_source: a.model,
      b_experiment: b.experiment, b_source: b.model,
    }));
  });
}

function setComparisonRunState(running) {
  const btn = document.getElementById('cmp-run-btn');
  if (!btn) return;
  btn.disabled    = running;
  btn.textContent = running ? '⏳ Comparing…' : '📊 Compare';
  btn.classList.toggle('running', running);
}

/**
 * Show a finished comparison  ← called by Python.
 * payload = { a:{label,stats}, b:{label,stats}, overlap, difference, units, plot_b64 }
 */
function setComparisonResult(payload) {
  const box  = document.getElementById('cmp-results');
  const body = document.getElementById('cmp-result-body');
  const img  = document.getElementById('cmp-plot');
  if (!box || !body) return;
  body.innerHTML = '';

  if (img) {
    if (payload.plot_b64) {
      img.src = 'data:image/png;base64,' + payload.plot_b64;
      img.style.display = 'block';
    } else {
      img.removeAttribute('src');
      img.style.display = 'none';
    }
  }

  const a = payload.a || {}, b = payload.b || {};
  const sa = a.stats || {}, sb = b.stats || {};
  body.appendChild(metricsTable(
    'Summary  [' + (payload.units || '—') + ']',
    ['', a.label || 'A', b.label || 'B'],
    [
      ['Mean',           sa.mean, sb.mean],
      ['Std dev',        sa.std,  sb.std],
      ['Min',            sa.min,  sb.min],
      ['Max',            sa.max,  sb.max],
      ['Trend / decade', sa.trend_per_decade, sb.trend_per_decade],
      ['Period',         sa.start + ' → ' + sa.end, sb.start + ' → ' + sb.end],
    ]));

  if (payload.difference !== null && payload.difference !== undefined) {
    const d = document.createElement('div');
    d.innerHTML = `Mean difference (A − B): <strong>${payload.difference}` +
                  `${payload.units ? ' ' + escapeHtml(payload.units) : ''}</strong>`;
    body.appendChild(d);
  }

  if (payload.overlap) {
    const o = payload.overlap;
    body.appendChild(metricsTable(
      `Overlap  ${o.start} → ${o.end}  (${o.n_days} steps)`,
      ['Metric', 'A vs B'],
      [['Bias', o.bias], ['RMSE', o.rmse], ['MAE', o.mae],
       ['Correlation', o.correlation]]));
  } else {
    const n = document.createElement('div');
    n.className = 'result-note';
    n.textContent = 'The two periods do not overlap — no paired metrics.';
    body.appendChild(n);
  }

  box.style.display = '';
}

/* ── Small result helpers ────────────────────────────────────────── */

/**
 * Build a table. Each row is [label, ...values]; a trailing 'abs' or 'high'
 * marks which of two numeric columns is the better one, so an improvement is
 * visible at a glance.
 */
function metricsTable(caption, headers, rows) {
  const table = document.createElement('table');
  table.className = 'result-table';

  const cap = document.createElement('caption');
  cap.textContent = caption;
  table.appendChild(cap);

  const thead = document.createElement('tr');
  headers.forEach(h => {
    const th = document.createElement('th');
    th.textContent = h;
    thead.appendChild(th);
  });
  table.appendChild(thead);

  rows.forEach(row => {
    const better = (row.length === headers.length + 1) ? row[row.length - 1] : null;
    const cells  = better ? row.slice(0, -1) : row;
    const tr = document.createElement('tr');
    cells.forEach((value, i) => {
      const td = document.createElement('td');
      td.textContent = (value === null || value === undefined || value === '')
        ? '—' : value;
      if (better && i === 2) {
        const prev = Number(cells[1]), cur = Number(value);
        // Identical values are neither an improvement nor a regression
        if (isFinite(prev) && isFinite(cur) && prev !== cur) {
          const improved = better === 'high' ? cur > prev
                                             : Math.abs(cur) < Math.abs(prev);
          td.classList.add(improved ? 'result-good' : 'result-bad');
        }
      }
      tr.appendChild(td);
    });
    table.appendChild(tr);
  });
  return table;
}

function escapeHtml(text) {
  const div = document.createElement('div');
  div.textContent = text;
  return div.innerHTML;
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
  ['tb-data-dir-display', 'basd-out-dir-display'].forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;
    el.textContent = path || 'default';
    el.title = path || 'Default from config';
  });
}

/* Opening defaults: an observational reference against a model, so both the
   BASD and Comparison tabs start on a combination that makes sense. */
const DATASET_DEFAULTS = {
  'basd-ref-dataset':    'mswx',
  'basd-target-dataset': 'cmip',
  'cmp-a-dataset':       'mswx',
  'cmp-b-dataset':       'cmip',
};

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
    const preferred = DATASET_DEFAULTS[id];
    if (cur && datasets.includes(cur))            sel.value = cur;
    else if (preferred && datasets.includes(preferred)) sel.value = preferred;
  });
}

/* ── CMIP6 model / experiment pickers ───────────────────────────── */

/* Every place a CMIP-driven dataset can be chosen. Python addresses these by
   slot name, so the same catalogue plumbing serves all four. */
const CMIP_SLOTS = {
  'download':    { dataset: 'tb-dataset',          section: 'cmip-options',
                   experiment: 'tb-experiment',    model: 'tb-model',
                   note: 'cmip-options-note' },
  'basd-target': { dataset: 'basd-target-dataset', section: 'cmip-options-basd-target',
                   experiment: 'basd-target-experiment', model: 'basd-target-model',
                   note: 'cmip-note-basd-target' },
  'cmp-a':       { dataset: 'cmp-a-dataset',       section: 'cmip-options-cmp-a',
                   experiment: 'cmp-a-experiment', model: 'cmp-a-model',
                   note: 'cmip-note-cmp-a' },
  'cmp-b':       { dataset: 'cmp-b-dataset',       section: 'cmip-options-cmp-b',
                   experiment: 'cmp-b-experiment', model: 'cmp-b-model',
                   note: 'cmip-note-cmp-b' },
};

/** True when *name* is a dataset configured by a CMIP6 model + experiment. */
function isCmipDataset(name) {
  return !!name && state.cmipDatasets.includes(String(name).toLowerCase());
}

/** Show one slot's model/experiment section only for CMIP-driven datasets. */
function updateCmipVisibility(slot, dataset) {
  const spec = CMIP_SLOTS[slot];
  if (!spec) return;
  if (dataset === undefined) {
    const ds = document.getElementById(spec.dataset);
    dataset = ds ? ds.value : '';
  }
  const section = document.getElementById(spec.section);
  if (section) section.style.display = isCmipDataset(dataset) ? '' : 'none';
}

/** Refresh visibility for every slot (after datasets are populated). */
function updateAllCmipVisibility() {
  Object.keys(CMIP_SLOTS).forEach(slot => updateCmipVisibility(slot));
}

/** Fill a <select> with *values*, keeping *selected* if present. */
function fillSelect(id, values, selected) {
  const sel = document.getElementById(id);
  if (!sel) return;
  sel.innerHTML = '';
  values.forEach(v => {
    const opt = document.createElement('option');
    opt.value = v; opt.textContent = v;
    sel.appendChild(opt);
  });
  if (selected && values.includes(selected)) sel.value = selected;
  sel.disabled = values.length === 0;
}

/**
 * Populate one slot's CMIP6 pickers  ← called by Python.
 * payload = { slot, experiments, experiment, models, model, loading, note }
 */
function setCmipOptions(payload) {
  const spec = CMIP_SLOTS[payload.slot || 'download'];
  if (!spec) return;
  const note = document.getElementById(spec.note);
  if (payload.loading) {
    [spec.experiment, spec.model].forEach(id => {
      const sel = document.getElementById(id);
      if (sel) sel.disabled = true;
    });
  } else {
    fillSelect(spec.experiment, payload.experiments || [], payload.experiment);
    fillSelect(spec.model,      payload.models      || [], payload.model);
  }
  if (note) {
    note.textContent = payload.note || '';
    note.classList.toggle('loading', !!payload.loading);
  }
}

/** Current model/experiment for a slot, or empty strings when hidden. */
function cmipSelection(slot) {
  const spec = CMIP_SLOTS[slot];
  const ds   = document.getElementById(spec.dataset);
  if (!ds || !isCmipDataset(ds.value)) return { experiment: '', model: '' };
  const exp = document.getElementById(spec.experiment);
  const mod = document.getElementById(spec.model);
  return { experiment: exp ? exp.value : '', model: mod ? mod.value : '' };
}

/** Wire the dataset / experiment / model selects of every slot. */
function initCmipSlots() {
  Object.entries(CMIP_SLOTS).forEach(([slot, spec]) => {
    const ds = document.getElementById(spec.dataset);
    if (ds && slot !== 'download') {
      // The download dataset select is handled in initToolbar (it also drives
      // the pipeline's dataset), the rest only steer their own pickers.
      ds.addEventListener('change', e => {
        updateCmipVisibility(slot, e.target.value);
        if (state.bridge) state.bridge.on_slot_dataset_changed(slot, e.target.value);
      });
    }
    const exp = document.getElementById(spec.experiment);
    if (exp) exp.addEventListener('change', e => {
      if (state.bridge) state.bridge.on_experiment_changed(slot, e.target.value);
    });
    const mod = document.getElementById(spec.model);
    if (mod) mod.addEventListener('change', e => {
      if (state.bridge) state.bridge.on_model_changed(slot, e.target.value);
    });
  });
}

function syncToolbar(payload) {
  if (payload.dataset) {
    const sel = document.getElementById('tb-dataset');
    if (sel && sel.querySelector(`option[value="${payload.dataset}"]`)) sel.value = payload.dataset;
    updateCmipVisibility('download', payload.dataset);
  }
  if (payload.start) document.getElementById('tb-start').value = payload.start;
  if (payload.end)   document.getElementById('tb-end').value   = payload.end;
  if (payload.data_dir !== undefined) setDataDir(payload.data_dir);
  if (payload.experiment) {
    const sel = document.getElementById('tb-experiment');
    if (sel && sel.querySelector(`option[value="${payload.experiment}"]`)) sel.value = payload.experiment;
  }
  if (payload.model) {
    const sel = document.getElementById('tb-model');
    if (sel && sel.querySelector(`option[value="${payload.model}"]`)) sel.value = payload.model;
  }
}

function onDashboardReady(payload) {
  if (Array.isArray(payload.cmip_datasets) && payload.cmip_datasets.length) {
    state.cmipDatasets = payload.cmip_datasets.map(d => String(d).toLowerCase());
  }
  populateDatasets(payload.datasets || []);
  syncToolbar(payload);
  updateAllCmipVisibility();
  // Confirm the actual displayed values back to Python so AppState stays in
  // sync even if this call races with an in-flight change message.
  const sel = document.getElementById('tb-dataset');
  if (sel && state.bridge) state.bridge.on_dataset_changed(sel.value);
  if (state.bridge) {
    Object.entries(CMIP_SLOTS).forEach(([slot, spec]) => {
      if (slot === 'download') return;
      const ds = document.getElementById(spec.dataset);
      if (ds) state.bridge.on_slot_dataset_changed(slot, ds.value);
    });
  }
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
      state.bridge.basd_run_state_changed.connect(setBasdRunState);
      state.bridge.comparison_run_state_changed.connect(setComparisonRunState);
      state.bridge.restarting.connect(setRestarting);

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
  initRestart();
  initBasemap();
  initControls();
  initOverlayControls();
  initMapEvents();
  initToolbar();
  initCmipSlots();
  initBasd();
  initComparison();
  initLog();
  initBridge();
})();

