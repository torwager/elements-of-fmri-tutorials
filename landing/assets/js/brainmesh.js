/* ==========================================================================
   Interactive triangulated brain mesh — home hero
   A lateral brain silhouette triangulated into vertices and edges.
   Moving the cursor near a vertex "fires" it: yellow impulses travel along
   mesh edges to connected vertices, which fire in turn (with refractory
   periods) so cascades ripple across the brain and die out naturally.
   Engine adapted from the neuron-network hero of torwager/mindfmricourse.
   ========================================================================== */
(function () {
  const canvas = document.getElementById('brainmesh');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // --- palette ------------------------------------------------------------
  const C = {
    edge: 'rgba(118,127,140,0.22)',
    edgeHot: '239,185,60',
    vert: '#c3c8cf',
    vertEdge: '#aab1ba',
    hot: '#f7c948',
    hotCore: '#fff3c4',
    membrane: '233,185,73',
  };

  // --- lateral brain silhouette (normalized, x∈[0,1], y∈[0,1], faces left)
  const OUTLINE = [
    [0.055, 0.42], [0.075, 0.33], [0.115, 0.245], [0.175, 0.17], [0.255, 0.11],
    [0.35, 0.065], [0.455, 0.04], [0.56, 0.04], [0.66, 0.06], [0.745, 0.10],
    [0.825, 0.16], [0.885, 0.235], [0.93, 0.32], [0.955, 0.41], [0.96, 0.50],
    [0.945, 0.575], [0.91, 0.635],                     /* occipital taper */
    [0.925, 0.69], [0.925, 0.76], [0.885, 0.825], [0.82, 0.86], [0.75, 0.855],
    [0.705, 0.82],                                     /* cerebellum */
    [0.685, 0.875], [0.665, 0.93], [0.635, 0.955], [0.615, 0.90], [0.62, 0.83],
    [0.60, 0.79],                                      /* brainstem */
    [0.54, 0.80], [0.46, 0.815], [0.375, 0.815], [0.29, 0.795], [0.215, 0.755],
    [0.165, 0.70], [0.14, 0.64],                       /* temporal underside */
    [0.10, 0.60], [0.065, 0.53], [0.05, 0.47],
  ];

  let W = 0, H = 0, dpr = 1;
  let nodes = [], edges = [], tris = [], impulses = [];
  let mouse = { x: -1e9, y: -1e9 };
  let running = true, last = performance.now(), ambientTimer = 1500;

  const rand = (a, b) => a + Math.random() * (b - a);

  function pointInPoly(x, y, poly) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
      if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }

  // --- Delaunay triangulation (Bowyer–Watson) -----------------------------
  function triangulate(pts) {
    const n = pts.length;
    // super-triangle
    let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9;
    for (const p of pts) { minx = Math.min(minx, p.x); miny = Math.min(miny, p.y); maxx = Math.max(maxx, p.x); maxy = Math.max(maxy, p.y); }
    const d = Math.max(maxx - minx, maxy - miny) * 10;
    const st = [{ x: minx - d, y: miny - d }, { x: maxx + d, y: miny - d }, { x: (minx + maxx) / 2, y: maxy + d }];
    const P = pts.concat(st);
    let triangles = [[n, n + 1, n + 2]];

    function circum(t) {
      const a = P[t[0]], b = P[t[1]], c = P[t[2]];
      const ax = a.x, ay = a.y, bx = b.x, by = b.y, cx = c.x, cy = c.y;
      const D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
      if (Math.abs(D) < 1e-12) return { x: 0, y: 0, r2: -1 };
      const ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / D;
      const uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / D;
      return { x: ux, y: uy, r2: (ax - ux) ** 2 + (ay - uy) ** 2 };
    }

    for (let i = 0; i < n; i++) {
      const p = P[i];
      const bad = [], polygon = [];
      for (const t of triangles) {
        const cc = circum(t);
        if (cc.r2 > 0 && (p.x - cc.x) ** 2 + (p.y - cc.y) ** 2 < cc.r2) bad.push(t);
      }
      for (const t of bad) {
        for (let e = 0; e < 3; e++) {
          const ed = [t[e], t[(e + 1) % 3]];
          let shared = false;
          for (const t2 of bad) {
            if (t2 === t) continue;
            for (let e2 = 0; e2 < 3; e2++) {
              const ed2 = [t2[e2], t2[(e2 + 1) % 3]];
              if ((ed[0] === ed2[0] && ed[1] === ed2[1]) || (ed[0] === ed2[1] && ed[1] === ed2[0])) { shared = true; break; }
            }
            if (shared) break;
          }
          if (!shared) polygon.push(ed);
        }
      }
      triangles = triangles.filter(t => !bad.includes(t));
      for (const ed of polygon) triangles.push([ed[0], ed[1], i]);
    }
    return triangles.filter(t => t[0] < n && t[1] < n && t[2] < n);
  }

  // --- build mesh ---------------------------------------------------------
  function build() {
    nodes = []; edges = []; tris = []; impulses = [];

    // fit the brain into the hero, biased right of center on wide screens
    const margin = 0.06 * Math.min(W, H);
    const availW = W > 900 ? W * 0.52 : W * 0.9;
    const availH = H * 0.78;
    const scale = Math.min(availW, availH * 1.28);
    const bw = scale, bh = scale / 1.28;
    const ox = W > 900 ? W * 0.44 + (W * 0.52 - bw) / 2 : (W - bw) / 2;
    const oy = (H - bh) / 2 - H * 0.02 + margin * 0;

    const poly = OUTLINE.map(([x, y]) => [ox + x * bw, oy + y * bh]);

    // vertices: outline points + jittered-grid interior points
    const spacing = Math.max(26, bw / (W < 700 ? 11 : 15));
    for (const [x, y] of poly) nodes.push({ x, y, r: rand(2.2, 3.4), light: 0, refr: 0, nbrs: [] });
    for (let gy = oy; gy < oy + bh; gy += spacing * 0.88) {
      for (let gx = ox; gx < ox + bw; gx += spacing) {
        const x = gx + ((Math.round((gy - oy) / (spacing * 0.88)) % 2) ? spacing / 2 : 0) + rand(-spacing * 0.26, spacing * 0.26);
        const y = gy + rand(-spacing * 0.24, spacing * 0.24);
        if (!pointInPoly(x, y, poly)) continue;
        // keep a little clearance from the outline so edge triangles stay clean
        let ok = true;
        for (const [px, py] of poly) { if ((px - x) ** 2 + (py - y) ** 2 < (spacing * 0.45) ** 2) { ok = false; break; } }
        if (!ok) continue;
        nodes.push({ x, y, r: rand(2.2, 3.8), light: 0, refr: 0, nbrs: [] });
      }
    }

    // Delaunay; drop triangles whose centroid leaves the silhouette
    const triangles = triangulate(nodes);
    const seen = new Set();
    for (const t of triangles) {
      const cx = (nodes[t[0]].x + nodes[t[1]].x + nodes[t[2]].x) / 3;
      const cy = (nodes[t[0]].y + nodes[t[1]].y + nodes[t[2]].y) / 3;
      if (!pointInPoly(cx, cy, poly)) continue;
      tris.push(t);
      for (let e = 0; e < 3; e++) {
        const i = t[e], j = t[(e + 1) % 3];
        const key = i < j ? i + '-' + j : j + '-' + i;
        if (seen.has(key)) continue;
        seen.add(key);
        const a = nodes[i], b = nodes[j];
        const ed = { a, b, ia: i, ib: j, hot: 0, len: Math.hypot(b.x - a.x, b.y - a.y) };
        edges.push(ed);
        a.nbrs.push({ e: ed, to: j }); b.nbrs.push({ e: ed, to: i });
      }
    }
  }

  function resize() {
    const rect = canvas.getBoundingClientRect();
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    W = rect.width; H = rect.height;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    build();
    if (!running) draw();
  }

  // --- firing -------------------------------------------------------------
  // Mesh vertices average ~6 neighbours (vs 2–3 in a sparse net), so
  // propagation probabilities are lower and fan-out is capped to keep
  // cascades local and organic.
  function fire(i, from, gen) {
    const n = nodes[i];
    if (n.refr > 0) return;
    n.light = 1; n.refr = 1100;
    if (reduceMotion) return;
    if (gen > 9) return;
    const cands = n.nbrs.filter(({ to }) => to !== from && nodes[to].refr <= 0);
    // shuffle, cap fan-out at 3
    for (let k = cands.length - 1; k > 0; k--) { const j = Math.floor(Math.random() * (k + 1)); [cands[k], cands[j]] = [cands[j], cands[k]]; }
    let sent = 0;
    for (const { e, to } of cands) {
      if (sent >= 3) break;
      const p = gen === 0 ? 0.95 : gen < 3 ? 0.62 : gen < 6 ? 0.42 : 0.25;
      if (Math.random() > p) continue;
      const forward = e.ia === i;
      impulses.push({ e, forward, t: 0, speed: rand(0.9, 1.4) * 260 / e.len, gen: gen + 1, from: i, to });
      sent++;
    }
  }

  function nearest(x, y, maxD) {
    let best = -1, bd = maxD * maxD;
    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i], d = (n.x - x) ** 2 + (n.y - y) ** 2;
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }

  // --- update -------------------------------------------------------------
  function update(dt) {
    const hi = nearest(mouse.x, mouse.y, 34);
    if (hi >= 0) fire(hi, -1, 0);

    if (!reduceMotion) {
      ambientTimer -= dt;
      if (ambientTimer <= 0 && impulses.length < 22) {
        fire(Math.floor(Math.random() * nodes.length), -1, 4);
        ambientTimer = rand(2200, 4800);
      }
    }

    for (const n of nodes) {
      if (n.refr > 0) n.refr -= dt;
      if (n.light > 0) n.light = Math.max(0, n.light - dt / 900);
    }
    for (const e of edges) if (e.hot > 0) e.hot = Math.max(0, e.hot - dt / 700);

    for (let k = impulses.length - 1; k >= 0; k--) {
      const im = impulses[k];
      im.t += im.speed * dt / 1000;
      im.e.hot = Math.max(im.e.hot, 0.9);
      if (im.t >= 1) { impulses.splice(k, 1); fire(im.to, im.from, im.gen); }
    }
  }

  // --- draw ---------------------------------------------------------------
  function lerpPoint(e, t) {
    const T = Math.min(1, Math.max(0, t));
    return { x: e.a.x + (e.b.x - e.a.x) * T, y: e.a.y + (e.b.y - e.a.y) * T };
  }

  function draw() {
    ctx.clearRect(0, 0, W, H);

    // membrane shimmer: faint amber fill on lit triangles
    for (const t of tris) {
      const avg = (nodes[t[0]].light + nodes[t[1]].light + nodes[t[2]].light) / 3;
      if (avg < 0.03) continue;
      ctx.beginPath();
      ctx.moveTo(nodes[t[0]].x, nodes[t[0]].y);
      ctx.lineTo(nodes[t[1]].x, nodes[t[1]].y);
      ctx.lineTo(nodes[t[2]].x, nodes[t[2]].y);
      ctx.closePath();
      ctx.fillStyle = `rgba(${C.membrane},${(avg * 0.08).toFixed(3)})`;
      ctx.fill();
    }

    // edges
    ctx.lineWidth = 1;
    for (const e of edges) {
      ctx.beginPath(); ctx.moveTo(e.a.x, e.a.y); ctx.lineTo(e.b.x, e.b.y);
      ctx.strokeStyle = C.edge; ctx.stroke();
      if (e.hot > 0.01) {
        ctx.strokeStyle = `rgba(${C.edgeHot},${(e.hot * 0.55).toFixed(3)})`;
        ctx.lineWidth = 1 + e.hot * 1.2; ctx.stroke(); ctx.lineWidth = 1;
      }
    }

    // vertices
    for (const n of nodes) {
      if (n.light > 0.02) {
        const g = ctx.createRadialGradient(n.x, n.y, 0, n.x, n.y, n.r + 20 * n.light);
        g.addColorStop(0, `rgba(255,214,90,${(0.55 * n.light).toFixed(3)})`);
        g.addColorStop(1, 'rgba(255,214,90,0)');
        ctx.fillStyle = g; ctx.beginPath(); ctx.arc(n.x, n.y, n.r + 20 * n.light, 0, Math.PI * 2); ctx.fill();
      }
      ctx.beginPath(); ctx.arc(n.x, n.y, n.r, 0, Math.PI * 2);
      if (n.light > 0.02) {
        const mix = n.light;
        ctx.fillStyle = `rgb(${Math.round(195 + 52 * mix)},${Math.round(200 + 1 * mix)},${Math.round(207 - 135 * mix)})`;
      } else ctx.fillStyle = C.vert;
      ctx.fill();
      ctx.strokeStyle = n.light > 0.02 ? C.hot : C.vertEdge; ctx.lineWidth = 1; ctx.stroke();
    }

    // impulses
    for (const im of impulses) {
      const T = im.forward ? im.t : 1 - im.t;
      const dir = im.forward ? 1 : -1;
      for (let s = 6; s >= 1; s--) {
        const p = lerpPoint(im.e, T - dir * s * 0.035);
        ctx.beginPath(); ctx.arc(p.x, p.y, 3.2 - s * 0.35, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(247,201,72,${(0.5 - s * 0.07).toFixed(3)})`; ctx.fill();
      }
      const p = lerpPoint(im.e, T);
      const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, 14);
      g.addColorStop(0, 'rgba(255,240,180,0.95)'); g.addColorStop(0.35, 'rgba(247,201,72,0.6)'); g.addColorStop(1, 'rgba(247,201,72,0)');
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(p.x, p.y, 14, 0, Math.PI * 2); ctx.fill();
      ctx.beginPath(); ctx.arc(p.x, p.y, 2.6, 0, Math.PI * 2); ctx.fillStyle = C.hotCore; ctx.fill();
    }
  }

  function loop(now) {
    if (!running) return;
    const dt = Math.min(50, now - last); last = now;
    update(dt); draw();
    requestAnimationFrame(loop);
  }

  // --- events -------------------------------------------------------------
  function toLocal(ev) {
    const r = canvas.getBoundingClientRect();
    return { x: ev.clientX - r.left, y: ev.clientY - r.top };
  }
  canvas.addEventListener('pointermove', ev => { mouse = toLocal(ev); }, { passive: true });
  canvas.addEventListener('pointerleave', () => { mouse = { x: -1e9, y: -1e9 }; });
  canvas.addEventListener('pointerdown', ev => {
    const p = toLocal(ev); const i = nearest(p.x, p.y, 44);
    if (i >= 0) { nodes[i].refr = 0; fire(i, -1, 0); }
  });
  const hero = canvas.parentElement;
  hero.addEventListener('pointermove', ev => { mouse = toLocal(ev); }, { passive: true });
  hero.addEventListener('pointerleave', () => { mouse = { x: -1e9, y: -1e9 }; });

  window.addEventListener('resize', () => { clearTimeout(window.__bmz); window.__bmz = setTimeout(resize, 120); });

  const io = new IntersectionObserver(entries => {
    const vis = entries[0].isIntersecting;
    if (vis && !running) { running = true; last = performance.now(); requestAnimationFrame(loop); }
    else if (!vis) running = false;
  }, { threshold: 0.02 });
  io.observe(canvas);

  resize();
  // Welcome cascade from somewhere near the brain's centre.
  setTimeout(() => { const i = nearest(W * 0.62, H * 0.45, 300); if (i >= 0) fire(i, -1, 0); }, 700);
  requestAnimationFrame(loop);
})();
