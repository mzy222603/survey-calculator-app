import { useState, useEffect, useCallback } from 'react';
import './App.css';

// 类型定义
interface Point { x: number; y: number; z?: number; name?: string; }
interface HistoryItem { expression: string; result: string; time: number; }
interface TraverseStation { angle: number; distance: number; }

// ==================== 测绘计算引擎 ====================
const Survey = {
  degToRad: (d: number) => d * Math.PI / 180,
  radToDeg: (r: number) => r * 180 / Math.PI,
  
  normalizeAz: (az: number) => { while(az<0)az+=360; while(az>=360)az-=360; return az; },
  
  dmsToD: (d: number, m: number, s: number) => {
    const sign = d >= 0 ? 1 : -1;
    return sign * (Math.abs(d) + m/60 + s/3600);
  },
  
  dToDms: (deg: number) => {
    const sign = deg >= 0 ? 1 : -1;
    deg = Math.abs(deg);
    const d = Math.floor(deg);
    const mf = (deg - d) * 60;
    const m = Math.floor(mf);
    const s = (mf - m) * 60;
    return { d: sign * d, m, s };
  },
  
  formatDms: (deg: number) => {
    const { d, m, s } = Survey.dToDms(deg);
    return `${d}°${m}'${s.toFixed(2)}"`;
  },
  
  // 坐标正算
  forward: (p: Point, az: number, dist: number): Point => ({
    x: p.x + dist * Math.cos(Survey.degToRad(az)),
    y: p.y + dist * Math.sin(Survey.degToRad(az))
  }),
  
  // 坐标反算
  inverse: (p1: Point, p2: Point) => {
    const dx = p2.x - p1.x, dy = p2.y - p1.y;
    const dist = Math.sqrt(dx*dx + dy*dy);
    let az = Survey.radToDeg(Math.atan2(dy, dx));
    return { azimuth: Survey.normalizeAz(az), distance: dist };
  },
  
  // 前方交会
  forwardIntersect: (pa: Point, pb: Point, angA: number, angB: number): Point => {
    const { azimuth: azAB } = Survey.inverse(pa, pb);
    const azAP = Survey.normalizeAz(azAB - angA);
    const azBP = Survey.normalizeAz(azAB + 180 + angB);
    const aR = Survey.degToRad(azAP), bR = Survey.degToRad(azBP);
    const denom = Math.sin(aR)*Math.cos(bR) - Math.cos(aR)*Math.sin(bR);
    const t = ((pb.y-pa.y)*Math.cos(bR) - (pb.x-pa.x)*Math.sin(bR)) / denom;
    return { x: pa.x + t*Math.cos(aR), y: pa.y + t*Math.sin(aR) };
  },
  
  // 后方交会
  resection: (pa: Point, pb: Point, pc: Point, alpha: number, beta: number): Point => {
    const a = Survey.degToRad(alpha), b = Survey.degToRad(beta);
    const cotA = 1/Math.tan(a), cotB = 1/Math.tan(b);
    const k1 = (pb.x-pa.x)*cotA - (pb.y-pa.y);
    const k2 = (pb.y-pa.y)*cotA + (pb.x-pa.x);
    const k3 = (pc.x-pb.x)*cotB - (pc.y-pb.y);
    const k4 = (pc.y-pb.y)*cotB + (pc.x-pb.x);
    const y = (k1-k3) / (k2-k4);
    const x = pa.x + (pa.y-y)*cotA;
    return { x, y };
  },
  
  // 多边形面积
  polyArea: (pts: Point[]) => {
    let a = 0;
    for(let i=0; i<pts.length; i++) {
      const j = (i+1) % pts.length;
      a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    return Math.abs(a) / 2;
  },
  
  // 闭合导线计算
  closedTraverse: (start: Point, startAz: number, stations: TraverseStation[]) => {
    const n = stations.length;
    const theory = (n) * 180; // 内角和理论值 = n*180 for closed
    const measured = stations.reduce((s, t) => s + t.angle, 0);
    const angClosure = measured - theory;
    const corr = -angClosure / n;
    
    const adjAngles = stations.map(s => s.angle + corr);
    const azimuths = [startAz];
    for(let i=0; i<n; i++) {
      azimuths.push(Survey.normalizeAz(azimuths[i] + adjAngles[i] - 180));
    }
    
    const dxList: number[] = [], dyList: number[] = [];
    for(let i=0; i<n; i++) {
      const az = Survey.degToRad(azimuths[i+1]);
      dxList.push(stations[i].distance * Math.cos(az));
      dyList.push(stations[i].distance * Math.sin(az));
    }
    
    const fx = dxList.reduce((a,b)=>a+b,0);
    const fy = dyList.reduce((a,b)=>a+b,0);
    const f = Math.sqrt(fx*fx + fy*fy);
    const totalLen = stations.reduce((s,t)=>s+t.distance,0);
    const relClosure = totalLen > 0 ? f/totalLen : 0;
    
    const points: Point[] = [start];
    let cx = start.x, cy = start.y;
    for(let i=0; i<n; i++) {
      const vx = -fx * stations[i].distance / totalLen;
      const vy = -fy * stations[i].distance / totalLen;
      cx += dxList[i] + vx;
      cy += dyList[i] + vy;
      points.push({ x: cx, y: cy, name: `P${i+1}` });
    }
    
    return { points, angClosure, fx, fy, f, relClosure, adjAngles };
  },
  
  // 附合导线
  attachedTraverse: (start: Point, end: Point, startAz: number, endAz: number, stations: TraverseStation[]) => {
    const n = stations.length;
    let calcAz = startAz;
    for(const s of stations) calcAz = Survey.normalizeAz(calcAz + s.angle - 180);
    
    const angClosure = Survey.normalizeAz(calcAz - endAz);
    const corr = -angClosure / n;
    const adjAngles = stations.map(s => s.angle + corr);
    
    const azimuths = [startAz];
    for(let i=0; i<n; i++) {
      azimuths.push(Survey.normalizeAz(azimuths[i] + adjAngles[i] - 180));
    }
    
    const dxList: number[] = [], dyList: number[] = [];
    for(let i=0; i<n; i++) {
      const az = Survey.degToRad(azimuths[i+1]);
      dxList.push(stations[i].distance * Math.cos(az));
      dyList.push(stations[i].distance * Math.sin(az));
    }
    
    const sumDx = dxList.reduce((a,b)=>a+b,0);
    const sumDy = dyList.reduce((a,b)=>a+b,0);
    const fx = sumDx - (end.x - start.x);
    const fy = sumDy - (end.y - start.y);
    const f = Math.sqrt(fx*fx + fy*fy);
    const totalLen = stations.reduce((s,t)=>s+t.distance,0);
    const relClosure = totalLen > 0 ? f/totalLen : 0;
    
    const points: Point[] = [start];
    let cx = start.x, cy = start.y, cumDist = 0;
    for(let i=0; i<n; i++) {
      cumDist += stations[i].distance;
      const vx = -fx * cumDist / totalLen;
      const vy = -fy * cumDist / totalLen;
      cx = start.x + dxList.slice(0,i+1).reduce((a,b)=>a+b,0) + vx;
      cy = start.y + dyList.slice(0,i+1).reduce((a,b)=>a+b,0) + vy;
      points.push({ x: cx, y: cy, name: `P${i+1}` });
    }
    
    return { points, angClosure, fx, fy, f, relClosure, adjAngles };
  },
  
  // 水准闭合路线平差
  levelClosed: (startH: number, diffs: number[], dists: number[]) => {
    const closure = diffs.reduce((a,b)=>a+b,0);
    const totalD = dists.reduce((a,b)=>a+b,0);
    const heights = [startH];
    for(let i=0; i<diffs.length; i++) {
      const v = -closure * dists[i] / totalD;
      heights.push(heights[i] + diffs[i] + v);
    }
    return { heights, closure };
  },
  
  // 水准附合路线平差
  levelAttached: (startH: number, endH: number, diffs: number[], dists: number[]) => {
    const measuredEnd = startH + diffs.reduce((a,b)=>a+b,0);
    const closure = measuredEnd - endH;
    const totalD = dists.reduce((a,b)=>a+b,0);
    const heights = [startH];
    for(let i=0; i<diffs.length; i++) {
      const v = -closure * dists[i] / totalD;
      heights.push(heights[i] + diffs[i] + v);
    }
    return { heights, closure };
  },
  
  // 高斯正算
  gaussForward: (lat: number, lon: number, L0?: number) => {
    const a = 6378137, f = 1/298.257222101;
    const e2 = 2*f - f*f;
    const zone = L0 ? Math.round((L0+3)/6) : Math.floor(lon/6)+1;
    const cm = L0 || zone*6-3;
    const B = Survey.degToRad(lat), l = Survey.degToRad(lon - cm);
    const e4=e2*e2, e6=e4*e2;
    const A0=1-e2/4-3*e4/64-5*e6/256;
    const X = a*(A0*B - 3/8*(e2+e4/4)*Math.sin(2*B) + 15/256*e4*Math.sin(4*B));
    const N = a/Math.sqrt(1-e2*Math.sin(B)*Math.sin(B));
    const t = Math.tan(B), t2=t*t;
    const eta2 = e2/(1-e2)*Math.cos(B)*Math.cos(B);
    const cB = Math.cos(B), l2=l*l;
    const x = X + N*t*cB*cB*l2/2*(1+(5-t2+9*eta2)*cB*cB*l2/12);
    let y = N*cB*l*(1+(1-t2+eta2)*cB*cB*l2/6) + 500000;
    return { x, y, zone, cm };
  },
  
  // 高斯反算
  gaussInverse: (x: number, y: number, cm: number) => {
    const a = 6378137, f = 1/298.257222101;
    const e2 = 2*f - f*f;
    y -= 500000;
    const e4=e2*e2, e6=e4*e2;
    const A0=1-e2/4-3*e4/64-5*e6/256;
    let Bf = x/(a*A0);
    for(let i=0; i<10; i++) {
      const FBf = a*(A0*Bf - 3/8*(e2+e4/4)*Math.sin(2*Bf) + 15/256*e4*Math.sin(4*Bf)) - x;
      const dF = a*A0*(1-e2*Math.sin(Bf)*Math.sin(Bf));
      Bf -= FBf/dF;
    }
    const Nf = a/Math.sqrt(1-e2*Math.sin(Bf)*Math.sin(Bf));
    const tf = Math.tan(Bf), tf2=tf*tf;
    const eta2f = e2/(1-e2)*Math.cos(Bf)*Math.cos(Bf);
    const B = Bf - tf*y*y/(2*Nf*Nf)*(1+eta2f);
    const l = y/(Nf*Math.cos(Bf))*(1-y*y/(6*Nf*Nf)*(1+2*tf2+eta2f));
    return { lat: Survey.radToDeg(B), lon: Survey.radToDeg(l)+cm };
  },
  
  // 四参数求解
  calc4Param: (src: Point[], tgt: Point[]) => {
    const n = src.length;
    let sXs=0,sYs=0,sXt=0,sYt=0,sXs2=0,sYs2=0,sXsXt=0,sYsYt=0,sXsYt=0,sYsXt=0;
    for(let i=0;i<n;i++) {
      sXs+=src[i].x; sYs+=src[i].y; sXt+=tgt[i].x; sYt+=tgt[i].y;
      sXs2+=src[i].x*src[i].x; sYs2+=src[i].y*src[i].y;
      sXsXt+=src[i].x*tgt[i].x; sYsYt+=src[i].y*tgt[i].y;
      sXsYt+=src[i].x*tgt[i].y; sYsXt+=src[i].y*tgt[i].x;
    }
    const A = sXs2+sYs2;
    const aa = (sXsXt+sYsYt)/A, bb = (sYsXt-sXsYt)/A;
    const dx = (sXt - aa*sXs + bb*sYs)/n;
    const dy = (sYt - aa*sYs - bb*sXs)/n;
    const scale = Math.sqrt(aa*aa+bb*bb);
    const rot = Math.atan2(bb,aa);
    return { dx, dy, scale, rotation: Survey.radToDeg(rot) };
  },
  
  // 四参数转换
  transform4: (pts: Point[], dx: number, dy: number, scale: number, rot: number) => {
    const r = Survey.degToRad(rot);
    const cosR = Math.cos(r), sinR = Math.sin(r);
    return pts.map(p => ({
      x: dx + scale*(p.x*cosR - p.y*sinR),
      y: dy + scale*(p.x*sinR + p.y*cosR)
    }));
  },
  
  // 圆曲线要素
  circularCurve: (R: number, alpha: number) => {
    const a = Survey.degToRad(Math.abs(alpha));
    return {
      T: R * Math.tan(a/2),
      L: R * a,
      E: R * (1/Math.cos(a/2) - 1),
      C: 2 * R * Math.sin(a/2)
    };
  },
  
  // 土方计算（断面法）
  earthwork: (areas: number[], dists: number[]) => {
    let vol = 0;
    for(let i=0; i<dists.length; i++) {
      vol += (areas[i] + areas[i+1]) / 2 * dists[i];
    }
    return vol;
  }
};

// ==================== 主应用 ====================
function App() {
  const [tab, setTab] = useState<'home'|'calc'|'survey'|'settings'>('home');
  const [display, setDisplay] = useState('0');
  const [expr, setExpr] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [mem, setMem] = useState(0);
  const [hasMem, setHasMem] = useState(false);
  const [angleUnit, setAngleUnit] = useState<'DEG'|'RAD'|'GRAD'>('DEG');
  const [precision, setPrecision] = useState(6);
  const [vibration, setVibration] = useState(true);
  const [theme, setTheme] = useState<'dark'|'light'>('dark');
  
  // 测绘计算状态
  const [surveyType, setSurveyType] = useState('forward');
  const [inputs, setInputs] = useState<{[k:string]:string}>({});
  const [result, setResult] = useState('');
  
  useEffect(() => {
    const saved = localStorage.getItem('survey_history');
    if(saved) setHistory(JSON.parse(saved));
  }, []);
  
  const vibrate = useCallback(() => {
    if(vibration && navigator.vibrate) navigator.vibrate(10);
  }, [vibration]);
  
  const fmt = (n: number) => n.toFixed(precision);
  
  const saveHistory = (e: string, r: string) => {
    const item = { expression: e, result: r, time: Date.now() };
    const h = [item, ...history].slice(0, 100);
    setHistory(h);
    localStorage.setItem('survey_history', JSON.stringify(h));
  };
  
  // 计算器函数
  const clear = () => { vibrate(); setDisplay('0'); setExpr(''); };
  const append = (v: string) => {
    vibrate();
    if(display === '0' && v !== '.') setDisplay(v);
    else if(display === 'Error') setDisplay(v);
    else setDisplay(display + v);
  };
  const back = () => { vibrate(); setDisplay(display.length > 1 ? display.slice(0,-1) : '0'); };
  const toggleSign = () => { vibrate(); if(display !== '0') setDisplay(display.startsWith('-') ? display.slice(1) : '-'+display); };
  
  const calc = () => {
    vibrate();
    try {
      let e = display
        .replace(/×/g,'*').replace(/÷/g,'/').replace(/π/g,`(${Math.PI})`).replace(/\^/g,'**')
        .replace(/√\(/g,'Math.sqrt(').replace(/∛\(/g,'Math.cbrt(')
        .replace(/sin\(/g,`Math.sin(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''}`)
        .replace(/cos\(/g,`Math.cos(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''}`)
        .replace(/tan\(/g,`Math.tan(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''}`)
        .replace(/asin\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.asin(`)
        .replace(/acos\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.acos(`)
        .replace(/atan\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.atan(`)
        .replace(/ln\(/g,'Math.log(').replace(/log\(/g,'Math.log10(')
        .replace(/abs\(/g,'Math.abs(').replace(/exp\(/g,'Math.exp(');
      const r = eval(e);
      const res = fmt(r);
      saveHistory(display, res);
      setExpr(display + ' =');
      setDisplay(res);
    } catch { setDisplay('Error'); }
  };
  
  const applyFn = (fn: string) => {
    vibrate();
    try {
      const v = parseFloat(display);
      let r: number;
      const toRad = angleUnit==='DEG' ? Math.PI/180 : angleUnit==='GRAD' ? Math.PI/200 : 1;
      const toDeg = angleUnit==='DEG' ? 180/Math.PI : angleUnit==='GRAD' ? 200/Math.PI : 1;
      switch(fn) {
        case 'sin': r = Math.sin(v*toRad); break;
        case 'cos': r = Math.cos(v*toRad); break;
        case 'tan': r = Math.tan(v*toRad); break;
        case 'asin': r = Math.asin(v)*toDeg; break;
        case 'acos': r = Math.acos(v)*toDeg; break;
        case 'atan': r = Math.atan(v)*toDeg; break;
        case 'sinh': r = Math.sinh(v); break;
        case 'cosh': r = Math.cosh(v); break;
        case 'tanh': r = Math.tanh(v); break;
        case 'ln': r = Math.log(v); break;
        case 'log': r = Math.log10(v); break;
        case '√': r = Math.sqrt(v); break;
        case '∛': r = Math.cbrt(v); break;
        case 'x²': r = v*v; break;
        case 'x³': r = v*v*v; break;
        case '1/x': r = 1/v; break;
        case 'n!': r = Array.from({length:Math.round(v)},(_, i)=>i+1).reduce((a,b)=>a*b,1); break;
        case 'abs': r = Math.abs(v); break;
        case '10ˣ': r = Math.pow(10,v); break;
        case 'eˣ': r = Math.exp(v); break;
        case '%': r = v/100; break;
        case 'π': r = Math.PI; break;
        case 'e': r = Math.E; break;
        case 'rand': r = Math.random(); break;
        default: return;
      }
      saveHistory(`${fn}(${display})`, fmt(r));
      setDisplay(fmt(r));
    } catch { setDisplay('Error'); }
  };
  
  const insertFn = (fn: string) => { vibrate(); setDisplay(display==='0'||display==='Error' ? fn+'(' : display+fn+'('); };
  
  // 度分秒转换
  const deg2dms = () => {
    vibrate();
    try {
      const deg = parseFloat(display);
      const sign = deg < 0 ? '-' : '';
      const abs = Math.abs(deg);
      const d = Math.floor(abs);
      const mf = (abs-d)*60;
      const m = Math.floor(mf);
      const s = (mf-m)*60;
      const r = `${sign}${d}°${m}'${s.toFixed(4)}"`;
      saveHistory(`D→DMS(${display})`, r);
      setDisplay(r);
    } catch { setDisplay('Error'); }
  };
  
  const dms2deg = () => {
    vibrate();
    try {
      const inp = display.trim();
      let deg = 0;
      if(inp.includes('°')) {
        const parts = inp.replace(/['"′″]/g,' ').replace('°',' ').trim().split(/\s+/);
        const sign = inp.startsWith('-') ? -1 : 1;
        const d = Math.abs(parseFloat(parts[0]))||0;
        const m = parseFloat(parts[1])||0;
        const s = parseFloat(parts[2])||0;
        deg = sign * (d + m/60 + s/3600);
      } else if(/^-?\d+\.\d{4,}$/.test(inp)) {
        const val = parseFloat(inp);
        const sign = val < 0 ? -1 : 1;
        const abs = Math.abs(val);
        const d = Math.floor(abs);
        const dec = abs - d;
        const mm = Math.floor(dec*100);
        const ss = (dec*100-mm)*100;
        deg = sign * (d + mm/60 + ss/3600);
      } else {
        deg = parseFloat(inp);
      }
      const r = fmt(deg);
      saveHistory(`DMS→D(${display})`, r);
      setDisplay(r);
    } catch { setDisplay('Error'); }
  };
  
  const deg2rad = () => { vibrate(); try { const r = fmt(parseFloat(display)*Math.PI/180); saveHistory(`D→R(${display})`,r); setDisplay(r); } catch { setDisplay('Error'); } };
  const rad2deg = () => { vibrate(); try { const r = fmt(parseFloat(display)*180/Math.PI); saveHistory(`R→D(${display})`,r); setDisplay(r); } catch { setDisplay('Error'); } };
  
  // 内存
  const mc = () => { vibrate(); setMem(0); setHasMem(false); };
  const mr = () => { vibrate(); setDisplay(String(mem)); };
  const mAdd = () => { vibrate(); setMem(mem + (parseFloat(display)||0)); setHasMem(true); };
  const mSub = () => { vibrate(); setMem(mem - (parseFloat(display)||0)); setHasMem(true); };
  
  // 测绘输入
  const inp = (k: string, v: string) => setInputs({...inputs, [k]: v});
  const getN = (k: string) => parseFloat(inputs[k]||'0');
  
  // 测绘计算
  const calcSurvey = () => {
    vibrate();
    try {
      let r = '';
      switch(surveyType) {
        case 'forward': {
          const p = Survey.forward({x:getN('x0'),y:getN('y0')}, getN('az'), getN('dist'));
          r = `【坐标正算结果】\nX = ${fmt(p.x)}\nY = ${fmt(p.y)}`;
          break;
        }
        case 'inverse': {
          const res = Survey.inverse({x:getN('x1'),y:getN('y1')}, {x:getN('x2'),y:getN('y2')});
          r = `【坐标反算结果】\n方位角 = ${fmt(res.azimuth)}° (${Survey.formatDms(res.azimuth)})\n距离 = ${fmt(res.distance)} m`;
          break;
        }
        case 'forward_intersect': {
          const p = Survey.forwardIntersect({x:getN('xa'),y:getN('ya')}, {x:getN('xb'),y:getN('yb')}, getN('angA'), getN('angB'));
          r = `【前方交会结果】\nXp = ${fmt(p.x)}\nYp = ${fmt(p.y)}`;
          break;
        }
        case 'resection': {
          const p = Survey.resection({x:getN('xa'),y:getN('ya')}, {x:getN('xb'),y:getN('yb')}, {x:getN('xc'),y:getN('yc')}, getN('alpha'), getN('beta'));
          r = `【后方交会结果】\nXp = ${fmt(p.x)}\nYp = ${fmt(p.y)}`;
          break;
        }
        case 'area': {
          const pts: Point[] = [];
          for(let i=1; i<=10; i++) {
            const x = inputs[`ax${i}`], y = inputs[`ay${i}`];
            if(x && y) pts.push({x:parseFloat(x), y:parseFloat(y)});
          }
          if(pts.length < 3) { r = '至少需要3个顶点'; break; }
          const area = Survey.polyArea(pts);
          r = `【面积计算结果】
顶点数: ${pts.length}
面积 = ${fmt(area)} m²
面积 = ${fmt(area/10000)} 公顷
面积 = ${fmt(area/666.67)} 亩`;
          break;
        }
        case 'closed_traverse': {
          const stations: TraverseStation[] = [];
          for(let i=1; i<=10; i++) {
            const ang = inputs[`tang${i}`], dist = inputs[`tdist${i}`];
            if(ang && dist) stations.push({angle:parseFloat(ang), distance:parseFloat(dist)});
          }
          if(stations.length < 3) { r = '至少需要3个测站'; break; }
          const tr = Survey.closedTraverse({x:getN('tx0'),y:getN('ty0')}, getN('taz0'), stations);
          r = '【闭合导线计算结果】\n\n角度闭合差: ' + fmt(tr.angClosure) + '" (' + Survey.formatDms(tr.angClosure/3600) + ')\nfx = ' + fmt(tr.fx) + ' m\nfy = ' + fmt(tr.fy) + ' m\n全长闭合差: ' + fmt(tr.f) + ' m\n相对闭合差: 1/' + Math.round(1/tr.relClosure) + '\n\n平差后坐标:\n' + tr.points.map((p,i) => (p.name||'起点') + ': X=' + fmt(p.x) + ', Y=' + fmt(p.y)).join('\n');
          break;
        }
        case 'attached_traverse': {
          const stations: TraverseStation[] = [];
          for(let i=1; i<=10; i++) {
            const ang = inputs[`atang${i}`], dist = inputs[`atdist${i}`];
            if(ang && dist) stations.push({angle:parseFloat(ang), distance:parseFloat(dist)});
          }
          if(stations.length < 1) { r = '至少需要1个测站'; break; }
          const tr = Survey.attachedTraverse(
            {x:getN('atx0'),y:getN('aty0')}, {x:getN('atxe'),y:getN('atye')},
            getN('ataz0'), getN('ataze'), stations
          );
          r = '【附合导线计算结果】\n\n角度闭合差: ' + fmt(tr.angClosure) + '"\nfx = ' + fmt(tr.fx) + ' m\nfy = ' + fmt(tr.fy) + ' m\n全长闭合差: ' + fmt(tr.f) + ' m\n相对闭合差: 1/' + Math.round(1/tr.relClosure) + '\n\n平差后坐标:\n' + tr.points.map((p,i) => (i===0?'起点':p.name) + ': X=' + fmt(p.x) + ', Y=' + fmt(p.y)).join('\n');
          break;
        }
        case 'level_closed': {
          const diffs: number[] = [], dists: number[] = [];
          for(let i=1; i<=10; i++) {
            const d = inputs[`ldiff${i}`], l = inputs[`ldist${i}`];
            if(d && l) { diffs.push(parseFloat(d)); dists.push(parseFloat(l)); }
          }
          if(diffs.length < 1) { r = '至少需要1段观测'; break; }
          const lv = Survey.levelClosed(getN('lh0'), diffs, dists);
          r = '【闭合水准路线平差】\n\n已知高程: ' + fmt(getN('lh0')) + ' m\n闭合差: ' + fmt(lv.closure*1000) + ' mm\n\n平差后高程:\n' + lv.heights.map((h,i) => '点' + i + ': H=' + fmt(h) + ' m').join('\n');
          break;
        }
        case 'level_attached': {
          const diffs: number[] = [], dists: number[] = [];
          for(let i=1; i<=10; i++) {
            const d = inputs[`aldiff${i}`], l = inputs[`aldist${i}`];
            if(d && l) { diffs.push(parseFloat(d)); dists.push(parseFloat(l)); }
          }
          if(diffs.length < 1) { r = '至少需要1段观测'; break; }
          const lv = Survey.levelAttached(getN('alh0'), getN('alhe'), diffs, dists);
          r = '【附合水准路线平差】\n\n起点高程: ' + fmt(getN('alh0')) + ' m\n终点高程: ' + fmt(getN('alhe')) + ' m\n闭合差: ' + fmt(lv.closure*1000) + ' mm\n\n平差后高程:\n' + lv.heights.map((h,i) => '点' + i + ': H=' + fmt(h) + ' m').join('\n');
          break;
        }
        case 'gauss_forward': {
          const g = Survey.gaussForward(getN('glat'), getN('glon'), getN('gcm')||undefined);
          r = '【高斯正算结果】\n\n输入:\n纬度 B = ' + fmt(getN('glat')) + '°\n经度 L = ' + fmt(getN('glon')) + '°\n\n输出:\nX = ' + fmt(g.x) + ' m\nY = ' + fmt(g.y) + ' m\n带号 = ' + g.zone + '\n中央子午线 = ' + g.cm + '°';
          break;
        }
        case 'gauss_inverse': {
          const g = Survey.gaussInverse(getN('gix'), getN('giy'), getN('gicm'));
          r = '【高斯反算结果】\n\n输入:\nX = ' + fmt(getN('gix')) + ' m\nY = ' + fmt(getN('giy')) + ' m\n中央子午线 = ' + getN('gicm') + '°\n\n输出:\n纬度 B = ' + fmt(g.lat) + '° (' + Survey.formatDms(g.lat) + ')\n经度 L = ' + fmt(g.lon) + '° (' + Survey.formatDms(g.lon) + ')';
          break;
        }
        case 'transform4': {
          const src: Point[] = [], tgt: Point[] = [];
          for(let i=1; i<=5; i++) {
            const sx=inputs[`t4sx${i}`], sy=inputs[`t4sy${i}`], tx=inputs[`t4tx${i}`], ty=inputs[`t4ty${i}`];
            if(sx&&sy&&tx&&ty) {
              src.push({x:parseFloat(sx),y:parseFloat(sy)});
              tgt.push({x:parseFloat(tx),y:parseFloat(ty)});
            }
          }
          if(src.length < 2) { r = '至少需要2个公共点'; break; }
          const p = Survey.calc4Param(src, tgt);
          r = '【四参数求解结果】\n\n公共点数: ' + src.length + '\n\n转换参数:\nΔX = ' + fmt(p.dx) + ' m\nΔY = ' + fmt(p.dy) + ' m\n尺度因子 K = ' + p.scale.toFixed(9) + '\n旋转角 θ = ' + fmt(p.rotation) + '° (' + Survey.formatDms(p.rotation) + ')';
          break;
        }
        case 'curve': {
          const c = Survey.circularCurve(getN('cR'), getN('cAlpha'));
          r = '【圆曲线要素计算】\n\n输入:\n半径 R = ' + fmt(getN('cR')) + ' m\n偏角 α = ' + fmt(getN('cAlpha')) + '°\n\n计算结果:\n切线长 T = ' + fmt(c.T) + ' m\n曲线长 L = ' + fmt(c.L) + ' m\n外矢距 E = ' + fmt(c.E) + ' m\n弦长 C = ' + fmt(c.C) + ' m';
          break;
        }
        case 'earthwork': {
          const areas: number[] = [], dists: number[] = [];
          for(let i=1; i<=10; i++) {
            const a = inputs[`ewa${i}`];
            if(a) areas.push(parseFloat(a));
            const d = inputs[`ewd${i}`];
            if(d) dists.push(parseFloat(d));
          }
          if(areas.length < 2 || dists.length < 1) { r = '至少需要2个断面和1个间距'; break; }
          const vol = Survey.earthwork(areas, dists);
          r = '【土方计算结果】\n\n断面数: ' + areas.length + '\n间距段数: ' + dists.length + '\n\n土方体积 = ' + fmt(vol) + ' m³';
          break;
        }
        default: r = '请选择计算类型';
      }
      setResult(r);
      saveHistory(surveyType, r.split('\n')[0]);
    } catch(e: any) {
      setResult(`计算错误: ${e.message || e}`);
    }
  };

  const surveyTypes = [
    { id: 'forward', name: '坐标正算', icon: '📍' },
    { id: 'inverse', name: '坐标反算', icon: '📏' },
    { id: 'forward_intersect', name: '前方交会', icon: '🔺' },
    { id: 'resection', name: '后方交会', icon: '🎯' },
    { id: 'area', name: '面积计算', icon: '⬛' },
    { id: 'closed_traverse', name: '闭合导线', icon: '🔄' },
    { id: 'attached_traverse', name: '附合导线', icon: '➡️' },
    { id: 'level_closed', name: '闭合水准', icon: '📊' },
    { id: 'level_attached', name: '附合水准', icon: '📈' },
    { id: 'gauss_forward', name: '高斯正算', icon: '🌍' },
    { id: 'gauss_inverse', name: '高斯反算', icon: '🗺️' },
    { id: 'transform4', name: '四参数转换', icon: '🔄' },
    { id: 'curve', name: '曲线计算', icon: '🛣️' },
    { id: 'earthwork', name: '土方计算', icon: '🏗️' },
  ];
  
  const InputField = ({label, k, placeholder}: {label: string; k: string; placeholder?: string}) => (
    <div className="input-row">
      <label>{label}</label>
      <input type="text" inputMode="decimal" value={inputs[k]||''} onChange={e=>inp(k,e.target.value)} placeholder={placeholder||'0'} />
    </div>
  );

  const renderSurveyInputs = () => {
    switch(surveyType) {
      case 'forward':
        return <><InputField label="起点X" k="x0"/><InputField label="起点Y" k="y0"/><InputField label="方位角(°)" k="az"/><InputField label="距离(m)" k="dist"/></>;
      case 'inverse':
        return <><InputField label="点1 X" k="x1"/><InputField label="点1 Y" k="y1"/><InputField label="点2 X" k="x2"/><InputField label="点2 Y" k="y2"/></>;
      case 'forward_intersect':
        return <><InputField label="A点X" k="xa"/><InputField label="A点Y" k="ya"/><InputField label="B点X" k="xb"/><InputField label="B点Y" k="yb"/><InputField label="∠PAB(°)" k="angA"/><InputField label="∠PBA(°)" k="angB"/></>;
      case 'resection':
        return <><InputField label="A点X" k="xa"/><InputField label="A点Y" k="ya"/><InputField label="B点X" k="xb"/><InputField label="B点Y" k="yb"/><InputField label="C点X" k="xc"/><InputField label="C点Y" k="yc"/><InputField label="∠APB(°)" k="alpha"/><InputField label="∠BPC(°)" k="beta"/></>;
      case 'area':
        return <div className="area-inputs">{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="point-row"><span>P{i}</span><input type="text" inputMode="decimal" value={inputs[`ax${i}`]||''} onChange={e=>inp(`ax${i}`,e.target.value)} placeholder="X"/><input type="text" inputMode="decimal" value={inputs[`ay${i}`]||''} onChange={e=>inp(`ay${i}`,e.target.value)} placeholder="Y"/></div>)}</div>;
      case 'closed_traverse':
        return <><InputField label="起点X" k="tx0"/><InputField label="起点Y" k="ty0"/><InputField label="起始方位角(°)" k="taz0"/><div className="table-header"><span>测站</span><span>水平角(°)</span><span>边长(m)</span></div>{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`tang${i}`]||''} onChange={e=>inp(`tang${i}`,e.target.value)} placeholder="角度"/><input type="text" inputMode="decimal" value={inputs[`tdist${i}`]||''} onChange={e=>inp(`tdist${i}`,e.target.value)} placeholder="边长"/></div>)}</>;
      case 'attached_traverse':
        return <><InputField label="起点X" k="atx0"/><InputField label="起点Y" k="aty0"/><InputField label="起始方位角(°)" k="ataz0"/><InputField label="终点X" k="atxe"/><InputField label="终点Y" k="atye"/><InputField label="终止方位角(°)" k="ataze"/><div className="table-header"><span>测站</span><span>水平角(°)</span><span>边长(m)</span></div>{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`atang${i}`]||''} onChange={e=>inp(`atang${i}`,e.target.value)} placeholder="角度"/><input type="text" inputMode="decimal" value={inputs[`atdist${i}`]||''} onChange={e=>inp(`atdist${i}`,e.target.value)} placeholder="边长"/></div>)}</>;
      case 'level_closed':
        return <><InputField label="起点高程(m)" k="lh0"/><div className="table-header"><span>段</span><span>高差(m)</span><span>距离(km)</span></div>{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`ldiff${i}`]||''} onChange={e=>inp(`ldiff${i}`,e.target.value)} placeholder="高差"/><input type="text" inputMode="decimal" value={inputs[`ldist${i}`]||''} onChange={e=>inp(`ldist${i}`,e.target.value)} placeholder="距离"/></div>)}</>;
      case 'level_attached':
        return <><InputField label="起点高程(m)" k="alh0"/><InputField label="终点高程(m)" k="alhe"/><div className="table-header"><span>段</span><span>高差(m)</span><span>距离(km)</span></div>{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`aldiff${i}`]||''} onChange={e=>inp(`aldiff${i}`,e.target.value)} placeholder="高差"/><input type="text" inputMode="decimal" value={inputs[`aldist${i}`]||''} onChange={e=>inp(`aldist${i}`,e.target.value)} placeholder="距离"/></div>)}</>;
      case 'gauss_forward':
        return <><InputField label="纬度B(°)" k="glat" placeholder="如 30.5"/><InputField label="经度L(°)" k="glon" placeholder="如 114.3"/><InputField label="中央子午线(°)" k="gcm" placeholder="自动计算"/></>;
      case 'gauss_inverse':
        return <><InputField label="X坐标(m)" k="gix"/><InputField label="Y坐标(m)" k="giy"/><InputField label="中央子午线(°)" k="gicm"/></>;
      case 'transform4':
        return <><div className="transform-header">公共点坐标（至少2个）</div><div className="table-header"><span>点</span><span>源X</span><span>源Y</span><span>目标X</span><span>目标Y</span></div>{[1,2,3,4,5].map(i=><div key={i} className="table-row-4"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`t4sx${i}`]||''} onChange={e=>inp(`t4sx${i}`,e.target.value)} placeholder="X"/><input type="text" inputMode="decimal" value={inputs[`t4sy${i}`]||''} onChange={e=>inp(`t4sy${i}`,e.target.value)} placeholder="Y"/><input type="text" inputMode="decimal" value={inputs[`t4tx${i}`]||''} onChange={e=>inp(`t4tx${i}`,e.target.value)} placeholder="X'"/><input type="text" inputMode="decimal" value={inputs[`t4ty${i}`]||''} onChange={e=>inp(`t4ty${i}`,e.target.value)} placeholder="Y'"/></div>)}</>;
      case 'curve':
        return <><InputField label="圆曲线半径R(m)" k="cR"/><InputField label="偏角α(°)" k="cAlpha"/></>;
      case 'earthwork':
        return <><div className="table-header"><span>断面</span><span>面积(m²)</span><span>间距(m)</span></div>{[1,2,3,4,5,6,7,8,9,10].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`ewa${i}`]||''} onChange={e=>inp(`ewa${i}`,e.target.value)} placeholder="面积"/><input type="text" inputMode="decimal" value={inputs[`ewd${i}`]||''} onChange={e=>inp(`ewd${i}`,e.target.value)} placeholder="间距"/></div>)}</>;
      default: return null;
    }
  };

  return (
    <div className={`app ${theme}`}>
      <div className="main-content">
        {/* 首页 */}
        {tab === 'home' && (
          <div className="home-page">
            <h1>测绘计算器Pro</h1>
            <p className="subtitle">专业测绘计算 · 向导式操作</p>
            <div className="quick-grid">
              <div className="quick-card" onClick={()=>{setTab('calc');}}>
                <span className="icon">🔢</span>
                <span>科学计算</span>
              </div>
              <div className="quick-card" onClick={()=>{setTab('survey');setSurveyType('forward');}}>
                <span className="icon">📍</span>
                <span>坐标正反算</span>
              </div>
              <div className="quick-card" onClick={()=>{setTab('survey');setSurveyType('closed_traverse');}}>
                <span className="icon">🔄</span>
                <span>导线计算</span>
              </div>
              <div className="quick-card" onClick={()=>{setTab('survey');setSurveyType('level_closed');}}>
                <span className="icon">📊</span>
                <span>水准平差</span>
              </div>
              <div className="quick-card" onClick={()=>{setTab('survey');setSurveyType('gauss_forward');}}>
                <span className="icon">🌍</span>
                <span>高斯投影</span>
              </div>
              <div className="quick-card" onClick={()=>{setTab('survey');setSurveyType('transform4');}}>
                <span className="icon">🔄</span>
                <span>坐标转换</span>
              </div>
            </div>
            <div className="history-section">
              <h3>历史记录</h3>
              {history.length === 0 ? <p className="empty">暂无记录</p> : (
                <div className="history-list">
                  {history.slice(0,10).map((h,i) => (
                    <div key={i} className="history-item">
                      <span className="expr">{h.expression}</span>
                      <span className="res">{h.result}</span>
                    </div>
                  ))}
                </div>
              )}
            </div>
          </div>
        )}

        {/* 科学计算器 */}
        {tab === 'calc' && (
          <div className="calc-page">
            <div className="display-area">
              <div className="expression">{expr}</div>
              <div className="display">{display}</div>
              <div className="status-bar">
                <span className={hasMem ? 'active' : ''}>M</span>
                <span>{angleUnit}</span>
              </div>
            </div>
            <div className="sci-panel">
              <div className="sci-row">
                <button onClick={()=>applyFn('sin')}>sin</button>
                <button onClick={()=>applyFn('cos')}>cos</button>
                <button onClick={()=>applyFn('tan')}>tan</button>
                <button onClick={()=>insertFn('sinh')}>sinh</button>
                <button onClick={()=>insertFn('cosh')}>cosh</button>
              </div>
              <div className="sci-row">
                <button onClick={()=>applyFn('asin')}>sin⁻¹</button>
                <button onClick={()=>applyFn('acos')}>cos⁻¹</button>
                <button onClick={()=>applyFn('atan')}>tan⁻¹</button>
                <button onClick={()=>applyFn('ln')}>ln</button>
                <button onClick={()=>applyFn('log')}>log</button>
              </div>
              <div className="sci-row">
                <button onClick={()=>applyFn('√')}>√</button>
                <button onClick={()=>applyFn('∛')}>∛</button>
                <button onClick={()=>applyFn('x²')}>x²</button>
                <button onClick={()=>applyFn('x³')}>x³</button>
                <button onClick={()=>append('^')}>^</button>
              </div>
              <div className="sci-row">
                <button onClick={()=>applyFn('n!')}>n!</button>
                <button onClick={()=>applyFn('1/x')}>1/x</button>
                <button onClick={()=>applyFn('abs')}>|x|</button>
                <button onClick={()=>applyFn('eˣ')}>eˣ</button>
                <button onClick={()=>applyFn('10ˣ')}>10ˣ</button>
              </div>
              <div className="sci-row">
                <button onClick={deg2dms}>D→DMS</button>
                <button onClick={dms2deg}>DMS→D</button>
                <button onClick={deg2rad}>D→R</button>
                <button onClick={rad2deg}>R→D</button>
                <button onClick={()=>append('°')}>°</button>
              </div>
              <div className="sci-row">
                <button onClick={mc}>MC</button>
                <button onClick={mr}>MR</button>
                <button onClick={mAdd}>M+</button>
                <button onClick={mSub}>M-</button>
                <button onClick={()=>setAngleUnit(angleUnit==='DEG'?'RAD':angleUnit==='RAD'?'GRAD':'DEG')}>{angleUnit}</button>
              </div>
            </div>
            <div className="num-panel">
              <div className="num-row">
                <button className="func" onClick={clear}>C</button>
                <button className="func" onClick={()=>setDisplay('0')}>CE</button>
                <button className="func" onClick={back}>⌫</button>
                <button className="op" onClick={()=>append('÷')}>÷</button>
              </div>
              <div className="num-row">
                <button onClick={()=>append('7')}>7</button>
                <button onClick={()=>append('8')}>8</button>
                <button onClick={()=>append('9')}>9</button>
                <button className="op" onClick={()=>append('×')}>×</button>
              </div>
              <div className="num-row">
                <button onClick={()=>append('4')}>4</button>
                <button onClick={()=>append('5')}>5</button>
                <button onClick={()=>append('6')}>6</button>
                <button className="op" onClick={()=>append('-')}>-</button>
              </div>
              <div className="num-row">
                <button onClick={()=>append('1')}>1</button>
                <button onClick={()=>append('2')}>2</button>
                <button onClick={()=>append('3')}>3</button>
                <button className="op" onClick={()=>append('+')}>+</button>
              </div>
              <div className="num-row">
                <button onClick={toggleSign}>±</button>
                <button onClick={()=>append('0')}>0</button>
                <button onClick={()=>append('.')}>.</button>
                <button className="eq" onClick={calc}>=</button>
              </div>
            </div>
          </div>
        )}

        {/* 测绘计算 */}
        {tab === 'survey' && (
          <div className="survey-page">
            <div className="type-selector">
              <select value={surveyType} onChange={e=>{setSurveyType(e.target.value);setInputs({});setResult('');}}>
                {surveyTypes.map(t => <option key={t.id} value={t.id}>{t.icon} {t.name}</option>)}
              </select>
            </div>
            <div className="survey-form">
              {renderSurveyInputs()}
            </div>
            <button className="calc-btn" onClick={calcSurvey}>计 算</button>
            {result && <div className="survey-result"><pre>{result}</pre></div>}
          </div>
        )}

        {/* 设置 */}
        {tab === 'settings' && (
          <div className="settings-page">
            <h2>设置</h2>
            <div className="setting-item">
              <span>角度单位</span>
              <select value={angleUnit} onChange={e=>setAngleUnit(e.target.value as any)}>
                <option value="DEG">度 (DEG)</option>
                <option value="RAD">弧度 (RAD)</option>
                <option value="GRAD">梯度 (GRAD)</option>
              </select>
            </div>
            <div className="setting-item">
              <span>计算精度</span>
              <select value={precision} onChange={e=>setPrecision(parseInt(e.target.value))}>
                <option value="4">4位小数</option>
                <option value="6">6位小数</option>
                <option value="8">8位小数</option>
                <option value="10">10位小数</option>
              </select>
            </div>
            <div className="setting-item">
              <span>震动反馈</span>
              <button className={`toggle ${vibration?'on':''}`} onClick={()=>setVibration(!vibration)}>{vibration?'开':'关'}</button>
            </div>
            <div className="setting-item">
              <span>主题</span>
              <button className={`toggle ${theme==='dark'?'on':''}`} onClick={()=>setTheme(theme==='dark'?'light':'dark')}>{theme==='dark'?'深色':'浅色'}</button>
            </div>
            <div className="setting-item">
              <span>清除历史</span>
              <button className="danger" onClick={()=>{setHistory([]);localStorage.removeItem('survey_history');}}>清除</button>
            </div>
            <div className="about">
              <p>测绘计算器Pro v2.0</p>
              <p>专业测绘计算解决方案</p>
            </div>
          </div>
        )}
      </div>

      {/* 底部导航 */}
      <nav className="bottom-nav">
        <button className={tab==='home'?'active':''} onClick={()=>setTab('home')}><span>🏠</span>首页</button>
        <button className={tab==='calc'?'active':''} onClick={()=>setTab('calc')}><span>🔢</span>计算器</button>
        <button className={tab==='survey'?'active':''} onClick={()=>setTab('survey')}><span>📐</span>测绘</button>
        <button className={tab==='settings'?'active':''} onClick={()=>setTab('settings')}><span>⚙️</span>设置</button>
      </nav>
    </div>
  );
}

export default App;
