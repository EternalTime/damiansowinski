(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const Dv = 0.1, V_MAX = 0.4;

  let N = 128;
  let f = 0.055, k = 0.062, ratio = 2.0;
  let stepsPerFrame = 8;
  let shading = false;
  let brushMode = 0;            // 0 = seed, 1 = cut
  let brushRadius = 0.02;       // in uv units
  let running = false, frameId = null;

  /* ── WebGL state ── */
  let canvas, gl;
  let simProg, brushProg, dispProg;
  let tex = [null, null], fbo = [null, null], cur = 0;
  let rampTex, sampNearest, sampLinear, vao;
  let brushQueue = [];          // {x0,y0,x1,y1,mode}
  let pointerDown = false, lastUV = null;

  /* ── Palette ramp (256×1 RGBA8) ── */
  const RAMP_SIZE = 256;
  const rampPix = new Uint8Array(RAMP_SIZE * 4);
  (function () {
    const stops = [
      _rgb('--near-black'),
      _rgb('--teal-dark'),
      _rgb('--teal-light'),
      _rgb('--cyan'),
      _rgb('--pink-light'),
      _rgb('--pink-dark'),
    ];
    for (let i = 0; i < RAMP_SIZE; i++) {
      const t = i / (RAMP_SIZE - 1);
      const ft = t * (stops.length - 1);
      const lo = Math.floor(ft), hi = Math.min(lo + 1, stops.length - 1);
      const a = ft - lo;
      rampPix[i*4]   = Math.round(stops[lo][0] + a*(stops[hi][0]-stops[lo][0]));
      rampPix[i*4+1] = Math.round(stops[lo][1] + a*(stops[hi][1]-stops[lo][1]));
      rampPix[i*4+2] = Math.round(stops[lo][2] + a*(stops[hi][2]-stops[lo][2]));
      rampPix[i*4+3] = 255;
    }
  })();

  /* ── Shaders ── */
  const VERT = `#version 300 es
  out vec2 v_uv;
  void main() {
    vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    v_uv = p;
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
  }`;

  const SIM_FRAG = `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform vec2  uTexel;
  uniform float uF, uK, uDu, uDv, uDt;
  in vec2 v_uv;
  out vec4 frag;
  void main() {
    vec2 s  = texture(uState, v_uv).rg;
    vec2 lap = texture(uState, v_uv + vec2( uTexel.x, 0.0)).rg
             + texture(uState, v_uv + vec2(-uTexel.x, 0.0)).rg
             + texture(uState, v_uv + vec2(0.0,  uTexel.y)).rg
             + texture(uState, v_uv + vec2(0.0, -uTexel.y)).rg
             - 4.0 * s;
    float uvv = s.x * s.y * s.y;
    float nu = s.x + uDt * (uDu * lap.x - uvv + uF * (1.0 - s.x));
    float nv = s.y + uDt * (uDv * lap.y + uvv - (uF + uK) * s.y);
    frag = vec4(clamp(nu, 0.0, 1.0), clamp(nv, 0.0, 1.0), 0.0, 1.0);
  }`;

  const BRUSH_FRAG = `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform vec2  uP0, uP1;
  uniform float uRadius;
  uniform int   uMode;
  uniform float uSeed;
  in vec2 v_uv;
  out vec4 frag;
  float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7)) + uSeed) * 43758.5453);
  }
  float distSeg(vec2 p, vec2 a, vec2 b) {
    vec2 ab = b - a;
    float t = clamp(dot(p - a, ab) / max(dot(ab, ab), 1e-12), 0.0, 1.0);
    return length(p - (a + t * ab));
  }
  void main() {
    vec2 s = texture(uState, v_uv).rg;
    float d = distSeg(v_uv, uP0, uP1);
    if (d < uRadius) {
      if (uMode == 0) {
        s.x = 0.8 + (hash(v_uv) - 0.5) * 0.05;
        s.y = 0.4 + (hash(v_uv + 7.0) - 0.5) * 0.05;
      } else {
        s.y *= smoothstep(0.0, 1.0, d / uRadius);
      }
    }
    frag = vec4(s, 0.0, 1.0);
  }`;

  const DISP_FRAG = `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform sampler2D uRamp;
  uniform vec2  uTexel;
  uniform float uVMax;
  uniform int   uShading;
  in vec2 v_uv;
  out vec4 frag;
  void main() {
    float v = texture(uState, v_uv).g;
    vec3 col = texture(uRamp, vec2(clamp(v / uVMax, 0.0, 1.0), 0.5)).rgb;
    if (uShading == 1) {
      float gx = texture(uState, v_uv + vec2(uTexel.x, 0.0)).g
               - texture(uState, v_uv - vec2(uTexel.x, 0.0)).g;
      float gy = texture(uState, v_uv + vec2(0.0, uTexel.y)).g
               - texture(uState, v_uv - vec2(0.0, uTexel.y)).g;
      vec3 n = normalize(vec3(-gx * 6.0, -gy * 6.0, 1.0));
      float light = 0.55 + 0.45 * dot(n, normalize(vec3(-0.4, -0.6, 0.7)));
      col *= light;
    }
    frag = vec4(col, 1.0);
  }`;

  function compile(type, src) {
    const sh = gl.createShader(type);
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
      throw new Error('shader: ' + gl.getShaderInfoLog(sh));
    }
    return sh;
  }
  function makeProg(fragSrc) {
    const p = gl.createProgram();
    gl.attachShader(p, compile(gl.VERTEX_SHADER, VERT));
    gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fragSrc));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error('link: ' + gl.getProgramInfoLog(p));
    }
    return p;
  }

  function makeStateTextures() {
    for (let i = 0; i < 2; i++) {
      if (tex[i]) { gl.deleteTexture(tex[i]); gl.deleteFramebuffer(fbo[i]); }
      tex[i] = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex[i]);
      gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RG16F, N, N);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      fbo[i] = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[i]);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex[i], 0);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }

  function initState() {
    if (!gl) return;
    for (let i = 0; i < 2; i++) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[i]);
      gl.viewport(0, 0, N, N);
      gl.clearColor(1.0, 0.0, 0.0, 1.0);   // U = 1, V = 0
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    brushQueue.length = 0;
  }

  function setResolution(newN) {
    N = newN;
    if (!gl) return;
    makeStateTextures();
    initState();
  }

  function initGL() {
    gl = canvas.getContext('webgl2', { antialias: false, alpha: false, preserveDrawingBuffer: false });
    if (!gl) return false;
    const ext = gl.getExtension('EXT_color_buffer_float')
             || gl.getExtension('EXT_color_buffer_half_float');
    if (!ext) { gl = null; return false; }

    simProg   = makeProg(SIM_FRAG);
    brushProg = makeProg(BRUSH_FRAG);
    dispProg  = makeProg(DISP_FRAG);

    vao = gl.createVertexArray();

    sampNearest = gl.createSampler();
    gl.samplerParameteri(sampNearest, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_WRAP_T, gl.REPEAT);
    sampLinear = gl.createSampler();
    gl.samplerParameteri(sampLinear, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_WRAP_T, gl.REPEAT);

    rampTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, rampTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, RAMP_SIZE, 1);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, RAMP_SIZE, 1, gl.RGBA, gl.UNSIGNED_BYTE, rampPix);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    makeStateTextures();
    initState();
    return true;
  }

  function drawFullscreen() {
    gl.bindVertexArray(vao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }

  function applyBrush(op) {
    const dst = 1 - cur;
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[dst]);
    gl.viewport(0, 0, N, N);
    gl.useProgram(brushProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampNearest);
    gl.uniform1i(gl.getUniformLocation(brushProg, 'uState'), 0);
    gl.uniform2f(gl.getUniformLocation(brushProg, 'uP0'), op.x0, op.y0);
    gl.uniform2f(gl.getUniformLocation(brushProg, 'uP1'), op.x1, op.y1);
    gl.uniform1f(gl.getUniformLocation(brushProg, 'uRadius'), brushRadius);
    gl.uniform1i(gl.getUniformLocation(brushProg, 'uMode'), op.mode);
    gl.uniform1f(gl.getUniformLocation(brushProg, 'uSeed'), Math.random() * 100.0);
    drawFullscreen();
    cur = dst;
  }

  function simStep(dt, Du) {
    const dst = 1 - cur;
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[dst]);
    gl.viewport(0, 0, N, N);
    gl.useProgram(simProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampNearest);
    gl.uniform1i(gl.getUniformLocation(simProg, 'uState'), 0);
    gl.uniform2f(gl.getUniformLocation(simProg, 'uTexel'), 1 / N, 1 / N);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uF'), f);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uK'), k);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uDu'), Du);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uDv'), Dv);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uDt'), dt);
    drawFullscreen();
    cur = dst;
  }

  function render() {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.useProgram(dispProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampLinear);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, rampTex);
    gl.bindSampler(1, null);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uState'), 0);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uRamp'), 1);
    gl.uniform2f(gl.getUniformLocation(dispProg, 'uTexel'), 1 / N, 1 / N);
    gl.uniform1f(gl.getUniformLocation(dispProg, 'uVMax'), V_MAX);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uShading'), shading ? 1 : 0);
    drawFullscreen();
  }

  function loop() {
    if (gl) {
      // flush brush strokes (works while paused too)
      for (let i = 0; i < brushQueue.length; i++) applyBrush(brushQueue[i]);
      brushQueue.length = 0;
      if (pointerDown && lastUV) {
        applyBrush({ x0: lastUV[0], y0: lastUV[1], x1: lastUV[0], y1: lastUV[1], mode: brushMode });
      }
      if (running) {
        const Du = ratio * Dv;
        const dt = Math.min(1.0, 0.24 / Math.max(Du, Dv));  // explicit-Euler stability guard
        for (let s = 0; s < stepsPerFrame; s++) simStep(dt, Du);
      }
      render();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'gs',
    title: 'Gray&ndash;Scott &mdash; Reaction&ndash;Diffusion',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="gsReset()">Reset</button><button class="applet-shell-header-btn" id="gs-pause-btn" onclick="gsTogglePause()">Pause</button>`,


    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Presets</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn gs-preset active" data-f="0.055" data-k="0.062" data-r="2.0">Coral</button>
          <button class="applet-shell-btn gs-preset"        data-f="0.078" data-k="0.061" data-r="2.0">Worms</button>
          <button class="applet-shell-btn gs-preset"        data-f="0.029" data-k="0.057" data-r="2.0">Maze</button>
          <button class="applet-shell-btn gs-preset"        data-f="0.039" data-k="0.058" data-r="2.0">Holes</button>
          <button class="applet-shell-btn gs-preset"        data-f="0.037" data-k="0.065" data-r="2.0">Mitosis</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Brush</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn gs-brush active" data-mode="0">Seed</button>
          <button class="applet-shell-btn gs-brush"        data-mode="1">Cut</button>
          <button class="applet-shell-btn" id="gs-shade">Shaded</button>
        </div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Fine</span>
          <input type="range" id="gs-brush-size" min="0.005" max="0.08" step="0.001" value="0.02">
          <span class="applet-shell-side">Broad</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Resolution</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn gs-res active" data-n="128">128</button>
          <button class="applet-shell-btn gs-res"        data-n="256">256</button>
          <button class="applet-shell-btn gs-res"        data-n="512">512</button>
          <button class="applet-shell-btn gs-res"        data-n="1024">1024</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Speed (steps/frame)</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="gs-speed" min="1" max="32" step="1" value="8">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Feed Rate f</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="gs-feed" min="0.01" max="0.10" step="0.0005" value="0.055">
          <span class="applet-shell-side">High</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Kill Rate k</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="gs-kill" min="0.03" max="0.07" step="0.0005" value="0.062">
          <span class="applet-shell-side">High</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Diffusion Ratio D<sub>u</sub>/D<sub>v</sub></div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">1&times;</span>
          <input type="range" id="gs-dratio" min="1.0" max="4.0" step="0.05" value="2.0">
          <span class="applet-shell-side">4&times;</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      canvas.style.cursor = 'crosshair';
      canvas.style.touchAction = 'none';

      if (!gl) {
        if (!initGL()) {
          const msg = document.createElement('div');
          msg.style.cssText = 'padding:24px;color:var(--text-bright);font-size:15px;';
          msg.textContent = 'WebGL2 with float render targets is not available in this browser.';
          canvas.parentNode.replaceChild(msg, canvas);
          return;
        }
      }

      function uvCoords(e) {
        const rect = canvas.getBoundingClientRect();
        return [
          (e.clientX - rect.left) / rect.width,
          1.0 - (e.clientY - rect.top) / rect.height,
        ];
      }
      canvas.addEventListener('pointerdown', function (e) {
        canvas.setPointerCapture(e.pointerId);
        pointerDown = true;
        lastUV = uvCoords(e);
        brushQueue.push({ x0: lastUV[0], y0: lastUV[1], x1: lastUV[0], y1: lastUV[1], mode: brushMode });
      });
      canvas.addEventListener('pointermove', function (e) {
        if (!pointerDown) return;
        const uv = uvCoords(e);
        brushQueue.push({ x0: lastUV[0], y0: lastUV[1], x1: uv[0], y1: uv[1], mode: brushMode });
        lastUV = uv;
      });
      const up = function () { pointerDown = false; lastUV = null; };
      canvas.addEventListener('pointerup', up);
      canvas.addEventListener('pointercancel', up);

      running = true;
      const pb = document.getElementById('gs-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      pointerDown = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('gs-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ canvas: c, S }) {
      canvas = c;   // shell has reset canvas.width/height; render() picks it up
    },
  });

  window.gsOpen  = () => shell.open();
  window.gsClose = () => shell.close();
  window.gsReset = function () { initState(); };
  window.gsTogglePause = function () {
    running = !running;
    const pb = document.getElementById('gs-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  document.getElementById('gs-feed').addEventListener('input', function() {
    f = parseFloat(this.value);
  });
  document.getElementById('gs-kill').addEventListener('input', function() {
    k = parseFloat(this.value);
  });
  document.getElementById('gs-dratio').addEventListener('input', function() {
    ratio = parseFloat(this.value);
  });
  document.getElementById('gs-speed').addEventListener('input', function() {
    stepsPerFrame = parseInt(this.value);
  });
  document.getElementById('gs-brush-size').addEventListener('input', function() {
    brushRadius = parseFloat(this.value);
  });
  document.getElementById('gs-shade').addEventListener('click', function() {
    shading = !shading;
    this.classList.toggle('active', shading);
  });

  document.querySelectorAll('.gs-brush').forEach(btn => {
    btn.addEventListener('click', function() {
      document.querySelectorAll('.gs-brush').forEach(b => b.classList.remove('active'));
      this.classList.add('active');
      brushMode = parseInt(this.dataset.mode);
    });
  });

  document.querySelectorAll('.gs-res').forEach(btn => {
    btn.addEventListener('click', function() {
      document.querySelectorAll('.gs-res').forEach(b => b.classList.remove('active'));
      this.classList.add('active');
      setResolution(parseInt(this.dataset.n));
    });
  });

  document.querySelectorAll('.gs-preset').forEach(btn => {
    btn.addEventListener('click', function() {
      document.querySelectorAll('.gs-preset').forEach(b => b.classList.remove('active'));
      this.classList.add('active');
      f = parseFloat(this.dataset.f);
      k = parseFloat(this.dataset.k);
      ratio = parseFloat(this.dataset.r);
      document.getElementById('gs-feed').value   = f;
      document.getElementById('gs-kill').value   = k;
      document.getElementById('gs-dratio').value = ratio;
    });
  });

})();
