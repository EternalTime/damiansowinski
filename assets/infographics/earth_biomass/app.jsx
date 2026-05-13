// The Biomass of Earth — isometric cube infographic
// Data: Bar-On, Phillips & Milo, PNAS 2018 (rounded Gt C)

const { useMemo, useState, useEffect, Fragment } = React;

// ─── Palette ─────────────────────────────────────────────────────────────────

const TOTAL_GTC = 545;

const PALETTES = {
  eukarya:  { top: '#88a872', right: '#6f9059', left: '#516e3c', edge: '#3d5429' },
  bacteria: { top: '#d18867', right: '#b56e4d', left: '#925233', edge: '#6b3a1f' },
  archaea:  { top: '#8e7099', right: '#755682', left: '#594064', edge: '#3f2a48' },

  plants:   { top: '#7d9d63', right: '#638250', left: '#496239', edge: '#354727' },
  fungi:    { top: '#d2ad5a', right: '#b89047', left: '#947133', edge: '#6d521f' },
  protists: { top: '#6c9c98', right: '#54817e', left: '#3e6562', edge: '#2b4845' },
  animals:  { top: '#b56f4f', right: '#985838', left: '#774127', edge: '#552c18' },
};

// ─── Data ─ Domain → Kingdom → Phylum/Environment ────────────────────────────
// Each node carries:
//   gtC      : its biomass in gigatons of carbon
//   palette  : colour key
//   shape    : optional [W, D, H] pile dimensions — auto when omitted
//   childUnit: Gt-per-cube used to scale its children's piles when expanded
//   children : sub-groups

const TREE = [
  { id: 'eukarya', label: 'Eukarya', gtC: 468, palette: 'eukarya',
    blurb: 'Cells with a nucleus — plants, animals, fungi, protists',
    shape: [12, 8, 5], childUnit: 1,
    children: [
      { id: 'plants', label: 'Plantae', gtC: 450, palette: 'plants',
        blurb: 'Trees, herbs and mosses — the bulk by far',
        shape: [10, 9, 5], childUnit: 10,
        children: [
          { id: 'angiosperms',  label: 'Flowering plants', gtC: 320, palette: 'plants',
            blurb: 'Angiosperms — most trees, all grasses & crops' },
          { id: 'gymnosperms',  label: 'Conifers',         gtC: 100, palette: 'plants',
            blurb: 'Gymnosperms — boreal & temperate forests' },
          { id: 'other-plants', label: 'Ferns & mosses',   gtC: 30,  palette: 'plants',
            blurb: 'Bryophytes, lycophytes & ferns' },
        ] },
      { id: 'fungi', label: 'Fungi', gtC: 12, palette: 'fungi',
        blurb: 'Underground mycelium dwarfs any visible mushroom',
        shape: [4, 3, 1], childUnit: 0.5,
        children: [
          { id: 'soil-fungi',  label: 'Soil fungi',       gtC: 10, palette: 'fungi',
            blurb: 'Mycorrhizal & saprotrophic networks',
            childUnit: 0.2,
            children: [
              { id: 'mycorrhizal',  label: 'Mycorrhizal', gtC: 7, palette: 'fungi',
                blurb: 'Plant-root symbionts trading nutrients for sugars' },
              { id: 'saprotrophic', label: 'Saprotrophic', gtC: 3, palette: 'fungi',
                blurb: 'Decomposers of leaf litter & dead wood' },
            ] },
          { id: 'other-fungi', label: 'Marine & other',   gtC: 2,  palette: 'fungi',
            blurb: 'Marine fungi, lichens, yeasts' },
        ] },
      { id: 'protists', label: 'Protists', gtC: 4, palette: 'protists',
        blurb: 'Single-celled eukaryotes — algae, amoebae, slime moulds',
        shape: [2, 2, 1], childUnit: 0.2,
        children: [
          { id: 'marine-protists', label: 'Marine protists', gtC: 2, palette: 'protists',
            blurb: 'Phytoplankton & marine protozoa' },
          { id: 'soil-protists',   label: 'Soil protists',   gtC: 2, palette: 'protists',
            blurb: 'Amoebae, ciliates & flagellates in soil' },
        ] },
      { id: 'animals', label: 'Animalia', gtC: 2, palette: 'animals',
        blurb: 'Multicellular movers — just 0.4 % of the biosphere',
        shape: [2, 1, 1], childUnit: 0.02,
        children: [
          { id: 'arthropods', label: 'Arthropods', gtC: 1.0, palette: 'animals',
            blurb: 'Insects, crustaceans, arachnids — half of all animal mass',
            childUnit: 0.01,
            children: [
              { id: 'crustaceans', label: 'Crustaceans', gtC: 0.6, palette: 'animals',
                blurb: 'Mostly oceanic copepods — drifting hordes' },
              { id: 'insects',     label: 'Insects',     gtC: 0.3, palette: 'animals',
                blurb: 'Termites & ants dominate by mass' },
              { id: 'arachnids',   label: 'Arachnids',   gtC: 0.1, palette: 'animals',
                blurb: 'Spiders, mites, ticks' },
            ] },
          { id: 'chordates',  label: 'Chordates',  gtC: 0.87, palette: 'animals',
            blurb: 'Animals with a spinal cord — fish, mammals, birds',
            childUnit: 0.01,
            children: [
              { id: 'fish',     label: 'Fish',     gtC: 0.7, palette: 'animals',
                blurb: 'Mesopelagic species dominate, mostly unseen' },
              { id: 'mammals',  label: 'Mammals',  gtC: 0.17, palette: 'animals',
                blurb: '~96 % of mammal mass is livestock or human — only ~4 % is wild',
                childUnit: 0.002,
                children: [
                  { id: 'livestock', label: 'Livestock', gtC: 0.1, palette: 'animals',
                    blurb: 'Cattle alone outweigh all wild mammals ~10×',
                    childUnit: 0.001,
                    children: [
                      { id: 'cattle',          label: 'Cattle',         gtC: 0.06,  palette: 'animals',
                        blurb: '~1 billion head of cattle worldwide' },
                      { id: 'other-livestock', label: 'Other livestock',gtC: 0.015, palette: 'animals',
                        blurb: 'Buffalo, horses, camels, other' },
                      { id: 'pigs',            label: 'Pigs',           gtC: 0.013, palette: 'animals' },
                      { id: 'sheep-goats',     label: 'Sheep & goats',  gtC: 0.012, palette: 'animals' },
                    ] },
                  { id: 'humans', label: 'Humans', gtC: 0.06, palette: 'animals',
                    blurb: '~8 billion of us, ~60 Mt C in total' },
                  { id: 'wild-mammals', label: 'Wild mammals', gtC: 0.007, palette: 'animals',
                    blurb: 'A fraction of pre-human levels (~0.02 Gt C, ~70 % lost)',
                    childUnit: 0.0005,
                    children: [
                      { id: 'marine-wild', label: 'Marine wild', gtC: 0.004, palette: 'animals',
                        blurb: 'Whales, dolphins, seals' },
                      { id: 'land-wild',   label: 'Land wild',   gtC: 0.003, palette: 'animals',
                        blurb: 'Deer, rodents, big cats — everything else' },
                    ] },
                ] },
              { id: 'birds',    label: 'Birds',    gtC: 0.007, palette: 'animals',
                blurb: 'Farmed poultry now outweigh all wild birds combined',
                childUnit: 0.0005,
                children: [
                  { id: 'poultry',    label: 'Poultry',    gtC: 0.005, palette: 'animals',
                    blurb: 'Mostly farmed chickens — ~25 billion alive at any time' },
                  { id: 'wild-birds', label: 'Wild birds', gtC: 0.002, palette: 'animals' },
                ] },
              { id: 'reptiles', label: 'Reptiles & amphibians', gtC: 0.002, palette: 'animals',
                blurb: 'Combined biomass barely registers against fish' },
            ] },
          { id: 'annelids',   label: 'Annelids',   gtC: 0.2, palette: 'animals',
            blurb: 'Earthworms & other segmented worms' },
          { id: 'molluscs',   label: 'Molluscs',   gtC: 0.2, palette: 'animals',
            blurb: 'Snails, clams, cephalopods' },
          { id: 'cnidarians', label: 'Cnidarians', gtC: 0.1, palette: 'animals',
            blurb: 'Corals, jellyfish, anemones' },
          { id: 'nematodes',  label: 'Nematodes',  gtC: 0.02, palette: 'animals',
            blurb: 'Tiny roundworms — outnumber every other animal' },
        ] },
    ] },
  { id: 'bacteria', label: 'Bacteria', gtC: 70, palette: 'bacteria',
    blurb: "Prokaryotes — Earth's most ancient and abundant cells",
    shape: [7, 5, 2], childUnit: 1,
    children: [
      { id: 'deep-bacteria',  label: 'Deep subsurface',  gtC: 60,  palette: 'bacteria',
        blurb: 'Kilometres below ground in rock & sediment',
        childUnit: 1,
        children: [
          { id: 'continental-deep', label: 'Continental crust', gtC: 30, palette: 'bacteria',
            blurb: 'Cells trapped in rock pores — slow metabolisms' },
          { id: 'marine-sediments', label: 'Marine sediments',  gtC: 30, palette: 'bacteria',
            blurb: 'Seafloor mud, kilometres of buried microbes' },
        ] },
      { id: 'soil-bacteria',  label: 'Soil',             gtC: 7,   palette: 'bacteria',
        blurb: 'Topsoil & the rhizosphere' },
      { id: 'plant-bacteria', label: 'Plant-associated', gtC: 1.6, palette: 'bacteria',
        blurb: 'Inside and on the surfaces of plants' },
      { id: 'marine-bacteria',label: 'Marine',           gtC: 1.3, palette: 'bacteria',
        blurb: 'Open-ocean and seafloor bacteria' },
    ] },
  { id: 'archaea', label: 'Archaea', gtC: 7, palette: 'archaea',
    blurb: 'Ancient single-celled life that thrives in extremes',
    shape: [7, 1, 1], childUnit: 0.2,
    children: [
      { id: 'deep-archaea',  label: 'Deep subsurface', gtC: 6,   palette: 'archaea',
        blurb: 'Kilometres underground, often anaerobic' },
      { id: 'soil-archaea',  label: 'Soil',            gtC: 0.5, palette: 'archaea',
        blurb: 'Ammonia oxidisers & methanogens' },
      { id: 'marine-archaea',label: 'Marine',          gtC: 0.3, palette: 'archaea',
        blurb: 'Deep-sea & sediment archaea' },
    ] },
];

const TOP_UNIT = 1; // Gt per cube at the top level

function findNode(id, nodes = TREE) {
  for (const n of nodes) {
    if (n.id === id) return n;
    if (n.children) {
      const f = findNode(id, n.children);
      if (f) return f;
    }
  }
  return null;
}

// ─── Isometric cube geometry ─────────────────────────────────────────────────

const ISO_C = 0.8660254;
const ISO_S = 0.5;

function proj(x, y, z, S) {
  return { px: (x - y) * ISO_C * S, py: (x + y) * ISO_S * S - z * S };
}

function Cube({ x, y, z, S, palette }) {
  const v = (dx, dy, dz) => proj(x + dx, y + dy, z + dz, S);
  const a = v(0,0,1), b = v(1,0,1), c = v(1,1,1), d = v(0,1,1);
  const e = v(1,0,0), f = v(1,1,0), g = v(0,1,0);
  const pt = p => `${p.px.toFixed(2)},${p.py.toFixed(2)}`;
  return (
    <g>
      <polygon points={`${pt(g)} ${pt(f)} ${pt(c)} ${pt(d)}`} fill={palette.left}  stroke={palette.edge} strokeWidth="0.4" strokeLinejoin="miter" />
      <polygon points={`${pt(e)} ${pt(f)} ${pt(c)} ${pt(b)}`} fill={palette.right} stroke={palette.edge} strokeWidth="0.4" strokeLinejoin="miter" />
      <polygon points={`${pt(a)} ${pt(b)} ${pt(c)} ${pt(d)}`} fill={palette.top}   stroke={palette.edge} strokeWidth="0.4" strokeLinejoin="miter" />
    </g>
  );
}

// Front-of-pile is high y, high x. Iterating z then y then x (all 0→max) means the
// last cells filled are the front-right of the topmost layer, so any partial top
// layer leaves its missing cubes facing the viewer.
function buildPileCoords(shape, fill) {
  const [W, D, H] = shape;
  const out = [];
  let i = 0;
  for (let z = 0; z < H; z++)
    for (let y = 0; y < D; y++)
      for (let x = 0; x < W; x++) {
        if (i++ < fill) out.push({ x, y, z });
      }
  return out;
}

function pileBounds(shape, S) {
  const [W, D, H] = shape;
  const corners = [
    proj(0, 0, 0, S), proj(W, 0, 0, S),
    proj(0, D, 0, S), proj(W, D, 0, S),
    proj(0, 0, H, S), proj(W, 0, H, S),
    proj(0, D, H, S), proj(W, D, H, S),
  ];
  const xs = corners.map(c => c.px), ys = corners.map(c => c.py);
  return {
    minX: Math.min(...xs) - 1, maxX: Math.max(...xs) + 1,
    minY: Math.min(...ys) - 1, maxY: Math.max(...ys) + 1,
  };
}

// Auto-pick a pile shape from a cube count when no manual shape is given.
function autoShape(count) {
  if (count <= 0) return [1, 1, 1];
  if (count <= 7) return [count, 1, 1];
  if (count <= 14) {
    const cols = Math.min(count, 5);
    const rows = Math.ceil(count / cols);
    return [cols, rows, 1];
  }
  let h = 1;
  if (count > 240) h = 5;
  else if (count > 120) h = 4;
  else if (count > 60) h = 3;
  else if (count > 25) h = 2;
  const perLayer = Math.ceil(count / h);
  const cols = Math.ceil(Math.sqrt(perLayer * 1.25));
  const rows = Math.ceil(perLayer / cols);
  return [cols, rows, h];
}

// ─── Pile component ──────────────────────────────────────────────────────────

function Pile({ node, unit, S, clickable, onClick, dim, asExpanded }) {
  const cubeCount = Math.max(1, Math.round(node.gtC / unit));
  const shape = node.shape || autoShape(cubeCount);
  const cubes = useMemo(() => buildPileCoords(shape, cubeCount), [shape.join(','), cubeCount]);
  const sorted = useMemo(
    () => [...cubes].sort((a, b) => (a.x + a.y + a.z) - (b.x + b.y + b.z) || a.z - b.z),
    [cubes]
  );
  const palette = PALETTES[node.palette];
  const b = pileBounds(shape, S);
  const w = b.maxX - b.minX, h = b.maxY - b.minY;
  const pct = (node.gtC / TOTAL_GTC) * 100;
  const pctStr = formatPct(pct);
  const val = formatValueDisplay(node.gtC);

  return (
    <div
      className={`pile-col ${dim ? 'dim' : ''} ${clickable ? 'clickable' : ''} ${asExpanded ? 'open' : ''}`}
      onClick={clickable ? (e) => { e.stopPropagation(); onClick(); } : undefined}
      style={{ '--pal-top': palette.top }}
    >
      <div className="pile-svg">
        <svg viewBox={`${b.minX} ${b.minY} ${w} ${h}`} width={w} height={h}
             style={{ display: 'block', overflow: 'visible' }}>
          <ellipse cx={(b.minX + b.maxX) / 2} cy={b.maxY - 2}
                   rx={(w / 2) * 0.92} ry={S * 0.35}
                   fill="rgba(40, 28, 14, 0.10)" />
          {sorted.map((cu, i) => <Cube key={i} {...cu} S={S} palette={palette} />)}
        </svg>
      </div>
      <div className="pile-cap">
        <div className="cap-line">
          <span className="swatch" style={{ background: palette.top }} />
          <span className="cap-label">{node.label}</span>
          {clickable && <span className="cap-plus" aria-hidden="true">+</span>}
        </div>
        <div className="cap-num">
          <span className="big">{val.num}</span>
          <span className="unit"> {val.unit}</span>
          <span className="pct"> · {pctStr}</span>
        </div>
        {node.blurb && <div className="cap-blurb">{node.blurb}</div>}
      </div>
    </div>
  );
}

function formatGtC(v) {
  if (v >= 10) return Math.round(v).toString();
  if (v >= 1) return v.toFixed(1).replace(/\.0$/, '');
  if (v >= 0.1) return v.toFixed(2).replace(/0$/, '');
  if (v >= 0.01) return v.toFixed(3).replace(/0+$/, '');
  if (v >= 0.001) return v.toFixed(4).replace(/0+$/, '');
  return v.toExponential(1);
}

// For pile captions: keep Gt for >= 0.1 Gt, switch to Mt / kt below so small
// numbers stay punchy instead of dissolving into decimals.
function formatValueDisplay(v) {
  if (v >= 0.1) return { num: formatGtC(v), unit: 'Gt\u00a0C' };
  if (v >= 0.001) return { num: `${parseFloat((v * 1000).toFixed(2))}`, unit: 'Mt\u00a0C' };
  return { num: `${parseFloat((v * 1e6).toFixed(2))}`, unit: 'kt\u00a0C' };
}

function formatPct(pct) {
  if (pct >= 1) return `${Math.round(pct)}%`;
  if (pct >= 0.1) return `${pct.toFixed(1)}%`;
  if (pct >= 0.01) return `${pct.toFixed(2)}%`;
  if (pct >= 0.001) return `${pct.toFixed(3)}%`;
  return '<0.001%';
}

// Cube scale label — swap to Mt/kt for legibility at deep levels.
// Returns { num, unit } so the caller can style them separately.
function formatCubeUnit(u) {
  if (u >= 0.1) {
    const num = u >= 1 ? `${u}` : `${u.toFixed(2).replace(/0+$/, '').replace(/\.$/, '')}`;
    return { num, unit: 'Gt\u00a0C' };
  }
  if (u >= 0.001) {
    const mt = parseFloat((u * 1000).toFixed(2));
    return { num: `${mt}`, unit: 'Mt\u00a0C' };
  }
  const kt = parseFloat((u * 1e6).toFixed(2));
  return { num: `${kt}`, unit: 'kt\u00a0C' };
}

// ─── Recursive slot ──────────────────────────────────────────────────────────
// A "slot" is one branch in the layout. If the current path expands this node,
// the slot renders that node's children inline; otherwise it renders a Pile.

function Slot({ node, unit, S, path, depth, onSelect, dim }) {
  const isOpen = path[depth] === node.id && node.children && node.children.length;
  const anySiblingChosen = path[depth] != null;
  const sibDim = anySiblingChosen && !isOpen;
  const effDim = dim || sibDim;

  if (isOpen) {
    const childUnit = node.childUnit || unit;
    const expVal = formatValueDisplay(node.gtC);
    return (
      <div className="exp" data-depth={depth}>
        <div className="exp-row">
          {node.children.map(child => (
            <Slot key={child.id} node={child} unit={childUnit} S={S}
                  path={path} depth={depth + 1} onSelect={onSelect} dim={false} />
          ))}
        </div>
        <button
          className="exp-tag"
          onClick={(e) => { e.stopPropagation(); onSelect(path.slice(0, depth)); }}
          style={{ '--pal-top': PALETTES[node.palette].top }}
        >
          <span className="exp-tag-bracket">⌐</span>
          <span className="exp-tag-label">{node.label}</span>
          <span className="exp-tag-meta">{expVal.num}&nbsp;{expVal.unit} · 1▢ = {formatUnit(childUnit)}</span>
          <span className="exp-tag-x">close</span>
        </button>
      </div>
    );
  }

  return (
    <Pile node={node} unit={unit} S={S}
          clickable={!!(node.children && node.children.length)}
          onClick={() => onSelect([...path.slice(0, depth), node.id])}
          dim={effDim} />
  );
}

function formatUnit(u) {
  const { num, unit } = formatCubeUnit(u);
  return `${num} ${unit.replace(/\u00a0/g, ' ')}`;
}

// ─── Breadcrumb ──────────────────────────────────────────────────────────────

function Breadcrumb({ path, onSelect }) {
  const crumbs = [{ label: 'All life', id: null }];
  for (let i = 0; i < path.length; i++) {
    const n = findNode(path[i]);
    crumbs.push({ label: n ? n.label : path[i], id: path[i], pal: n && PALETTES[n.palette] });
  }
  return (
    <nav className="bc" aria-label="Selection">
      {crumbs.map((c, i) => (
        <span className="bc-group" key={i}>
          {i > 0 && <span className="bc-sep">›</span>}
          <button
            className={`bc-item ${i === crumbs.length - 1 ? 'on' : ''}`}
            onClick={() => onSelect(path.slice(0, i))}
            disabled={i === crumbs.length - 1}
          >
            {c.pal && <span className="bc-sw" style={{ background: c.pal.top }} />}
            {c.label}
          </button>
        </span>
      ))}
      {path.length > 0 && (
        <button className="bc-reset" onClick={() => onSelect([])}>Reset view</button>
      )}
    </nav>
  );
}

// ─── Tweaks-driven shell ─────────────────────────────────────────────────────

const TWEAK_DEFAULTS = /*EDITMODE-BEGIN*/{
  "cubeSize": 15,
  "showGround": true
}/*EDITMODE-END*/;

function App() {
  const [t, setTweak] = useTweaks(TWEAK_DEFAULTS);
  const [path, setPath] = useState([]);

  // ESC backs out one level
  useEffect(() => {
    const onKey = (e) => {
      if (e.key === 'Escape' && path.length) setPath(p => p.slice(0, -1));
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [path.length]);

  const S = t.cubeSize;

  // Compute a helpful caption for the current level's unit
  const currentUnit = useMemo(() => {
    let u = TOP_UNIT;
    let nodes = TREE;
    for (const id of path) {
      const n = nodes.find(x => x.id === id);
      if (!n) break;
      u = n.childUnit || u;
      nodes = n.children || [];
    }
    return u;
  }, [path]);

  return (
    <div className="page">
      <header className="hd">
        <div className="hd-text">
          <div className="eyebrow">An infographic · biosphere</div>
          <h1>The biomass<br/>of Earth</h1>
        </div>
        <div className="hd-meta">
          <p className="sub">
            All <em>living</em> carbon on the planet — ~{TOTAL_GTC}&nbsp;Gt C in total.
            Each cube is a gigaton of carbon. Click a stack to break it down.
          </p>
          <div className="legend">
            <span className="leg-k">1 cube</span>
            <span className="leg-v">= {formatUnit(currentUnit)}</span>
          </div>
        </div>
      </header>

      <div className="bc-wrap">
        <Breadcrumb path={path} onSelect={setPath} />
      </div>

      <section className={`piles depth-${path.length}`}>
        {t.showGround && <div className="ground" />}
        <div className="piles-row">
          {TREE.map(node => (
            <Slot key={node.id} node={node} unit={TOP_UNIT} S={S}
                  path={path} depth={0} onSelect={setPath} dim={false} />
          ))}
        </div>
      </section>

      <footer className="ft">
        <div className="ft-meta">
          <div className="ft-row"><span className="ft-k">Source</span>
            <span className="ft-v">Bar-On, Phillips &amp; Milo, <em>PNAS</em> 2018 — figures rounded.</span></div>
          <div className="ft-row"><span className="ft-k">Not pictured</span>
            <span className="ft-v">Viruses ~0.2 Gt C · dead biomass &amp; detritus excluded.</span></div>
          <div className="ft-row"><span className="ft-k">Interaction</span>
            <span className="ft-v">Click any stack to break it down — three levels deep. Use the breadcrumb, the close tag, or <kbd>Esc</kbd> to step back.</span></div>
        </div>
        <div className="ft-note">
          The carbon in <span className="mono">one cube</span> is roughly what
          <span className="mono"> 40&nbsp;billion</span> mature trees lock away each year.
          Earth's plants weigh in <span className="mono">200×</span> more than every animal combined.
        </div>
      </footer>

      <TweaksPanel>
        <TweakSection label="Cubes" />
        <TweakSlider label="Cube size" value={t.cubeSize} min={10} max={22} step={1} unit="px"
                     onChange={(v) => setTweak('cubeSize', v)} />
        <TweakToggle label="Show ground line" value={t.showGround}
                     onChange={(v) => setTweak('showGround', v)} />
      </TweaksPanel>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById('root')).render(<App />);
