import { useState, useEffect, useCallback, useRef } from 'react';
import './App.css';

// 类型定义
interface Point { x: number; y: number; z?: number; name?: string; }
interface HistoryItem { 
  expression: string; 
  result: string; 
  time: number; 
  category: 'calc' | 'survey'; 
  surveyType?: string;
  inputs?: {[k:string]:string};
}
interface TraverseStation { angle: number; distance: number; }

// 主题配色
const Themes: {[k:string]:{name:string;bg:string;card:string;primary:string;text:string;border:string}} = {
  'dark': { name: '深色', bg: '#0d1117', card: '#161b22', primary: '#238636', text: '#e6edf3', border: '#30363d' },
  'light': { name: '浅色', bg: '#f6f8fa', card: '#ffffff', primary: '#1a7f37', text: '#24292f', border: '#d0d7de' },
  'blue': { name: '蓝色', bg: '#0a1628', card: '#0f2744', primary: '#1f6feb', text: '#c9d1d9', border: '#21436d' },
  'green': { name: '绿色', bg: '#0d1810', card: '#132318', primary: '#2ea043', text: '#c5d9c8', border: '#1f3d25' },
  'purple': { name: '紫色', bg: '#150d1e', card: '#1c1329', primary: '#8957e5', text: '#d2c9e0', border: '#3d2a5c' },
  'orange': { name: '橙色', bg: '#1a1008', card: '#241a10', primary: '#d29922', text: '#e0d4c0', border: '#4a3520' }
};

// ==================== 测绘计算引擎 ====================
// 椭球参数
const Ellipsoids: {[k:string]:{a:number;f:number;name:string}} = {
  'CGCS2000': { a: 6378137.0, f: 1/298.257222101, name: 'CGCS2000/国家大地坐标系' },
  'WGS84': { a: 6378137.0, f: 1/298.257223563, name: 'WGS84/GPS' },
  'BJ54': { a: 6378245.0, f: 1/298.3, name: '北京54' },
  'XIAN80': { a: 6378140.0, f: 1/298.257, name: '西安80' },
  'KRASOVSKY': { a: 6378245.0, f: 1/298.3, name: '克拉索夫斯基' },
  'GRS80': { a: 6378137.0, f: 1/298.257222101, name: 'GRS80' }
};

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
  },

  // ========== 全面坐标转换功能 ==========
  
  // 大地坐标(BLH) → 空间直角坐标(XYZ)
  blhToXyz: (B: number, L: number, H: number, ellipsoid: string = 'CGCS2000') => {
    const { a, f } = Ellipsoids[ellipsoid] || Ellipsoids['CGCS2000'];
    const e2 = 2*f - f*f;
    const Br = Survey.degToRad(B), Lr = Survey.degToRad(L);
    const sinB = Math.sin(Br), cosB = Math.cos(Br);
    const N = a / Math.sqrt(1 - e2 * sinB * sinB);
    return {
      X: (N + H) * cosB * Math.cos(Lr),
      Y: (N + H) * cosB * Math.sin(Lr),
      Z: (N * (1 - e2) + H) * sinB
    };
  },
  
  // 空间直角坐标(XYZ) → 大地坐标(BLH)
  xyzToBlh: (X: number, Y: number, Z: number, ellipsoid: string = 'CGCS2000') => {
    const { a, f } = Ellipsoids[ellipsoid] || Ellipsoids['CGCS2000'];
    const e2 = 2*f - f*f;
    const b = a * (1 - f);
    const ep2 = (a*a - b*b) / (b*b);
    const p = Math.sqrt(X*X + Y*Y);
    const theta = Math.atan2(Z * a, p * b);
    const L = Math.atan2(Y, X);
    const B = Math.atan2(Z + ep2 * b * Math.pow(Math.sin(theta), 3), p - e2 * a * Math.pow(Math.cos(theta), 3));
    const sinB = Math.sin(B);
    const N = a / Math.sqrt(1 - e2 * sinB * sinB);
    const H = p / Math.cos(B) - N;
    return { B: Survey.radToDeg(B), L: Survey.radToDeg(L), H };
  },
  
  // 七参数转换（布尔萨模型）
  // 参数: dx,dy,dz(平移m), rx,ry,rz(旋转角秒), m(尺度ppm)
  transform7Param: (X: number, Y: number, Z: number, dx: number, dy: number, dz: number, rx: number, ry: number, rz: number, m: number) => {
    // 角度从角秒转弧度
    const rxRad = rx * Math.PI / 648000; // 角秒转弧度
    const ryRad = ry * Math.PI / 648000;
    const rzRad = rz * Math.PI / 648000;
    const scale = 1 + m * 1e-6; // ppm转尺度因子
    
    // 布尔萨公式
    const Xn = dx + scale * (X - rzRad * Y + ryRad * Z);
    const Yn = dy + scale * (rzRad * X + Y - rxRad * Z);
    const Zn = dz + scale * (-ryRad * X + rxRad * Y + Z);
    return { X: Xn, Y: Yn, Z: Zn };
  },
  
  // 七参数求解（至少3个公共点）
  calc7Param: (src: {X:number;Y:number;Z:number}[], tgt: {X:number;Y:number;Z:number}[]) => {
    const n = src.length;
    if (n < 3) throw new Error('至少需要3个公共点');
    
    // 简化的最小二乘法求解
    let sumDx = 0, sumDy = 0, sumDz = 0;
    for (let i = 0; i < n; i++) {
      sumDx += tgt[i].X - src[i].X;
      sumDy += tgt[i].Y - src[i].Y;
      sumDz += tgt[i].Z - src[i].Z;
    }
    const dx = sumDx / n, dy = sumDy / n, dz = sumDz / n;
    
    // 简化计算旋转和尺度
    let sumScale = 0, count = 0;
    for (let i = 0; i < n; i++) {
      const srcLen = Math.sqrt(src[i].X*src[i].X + src[i].Y*src[i].Y + src[i].Z*src[i].Z);
      const tgtLen = Math.sqrt(tgt[i].X*tgt[i].X + tgt[i].Y*tgt[i].Y + tgt[i].Z*tgt[i].Z);
      if (srcLen > 0) { sumScale += tgtLen / srcLen; count++; }
    }
    const m = count > 0 ? (sumScale / count - 1) * 1e6 : 0;
    
    return { dx, dy, dz, rx: 0, ry: 0, rz: 0, m };
  },
  
  // 三参数转换（仅平移）
  transform3Param: (X: number, Y: number, Z: number, dx: number, dy: number, dz: number) => ({
    X: X + dx, Y: Y + dy, Z: Z + dz
  }),
  
  // 高斯投影（支持不同椭球和带宽）
  gaussProj: (B: number, L: number, zoneWidth: 3|6, L0?: number, ellipsoid: string = 'CGCS2000') => {
    const { a, f } = Ellipsoids[ellipsoid] || Ellipsoids['CGCS2000'];
    const e2 = 2*f - f*f;
    
    // 计算带号和中央子午线
    let zone: number, cm: number;
    if (L0 !== undefined) {
      cm = L0;
      zone = zoneWidth === 6 ? Math.floor((L0 + 6) / 6) : Math.floor((L0 + 1.5) / 3);
    } else {
      if (zoneWidth === 6) {
        zone = Math.floor(L / 6) + 1;
        cm = zone * 6 - 3;
      } else {
        zone = Math.floor((L + 1.5) / 3);
        cm = zone * 3;
      }
    }
    
    const Br = Survey.degToRad(B), l = Survey.degToRad(L - cm);
    const e4 = e2*e2, e6 = e4*e2;
    const A0 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
    const A2 = 3/8 * (e2 + e4/4 + 15*e6/128);
    const A4 = 15/256 * (e4 + 3*e6/4);
    const A6 = 35*e6/3072;
    const X0 = a * (A0*Br - A2*Math.sin(2*Br) + A4*Math.sin(4*Br) - A6*Math.sin(6*Br));
    
    const sinB = Math.sin(Br), cosB = Math.cos(Br);
    const N = a / Math.sqrt(1 - e2*sinB*sinB);
    const t = Math.tan(Br), t2 = t*t;
    const eta2 = e2/(1-e2) * cosB*cosB;
    const l2 = l*l;
    
    const x = X0 + N*t*cosB*cosB*l2/2 * (1 + (5-t2+9*eta2+4*eta2*eta2)*cosB*cosB*l2/12 + (61-58*t2+t2*t2)*Math.pow(cosB*l,4)/360);
    const y = N*cosB*l * (1 + (1-t2+eta2)*cosB*cosB*l2/6 + (5-18*t2+t2*t2+14*eta2-58*eta2*t2)*Math.pow(cosB*l,4)/120);
    
    return { x, y: y + 500000, zone, cm, zoneWidth };
  },
  
  // 高斯反算（支持不同椭球）
  gaussInvProj: (x: number, y: number, cm: number, ellipsoid: string = 'CGCS2000') => {
    const { a, f } = Ellipsoids[ellipsoid] || Ellipsoids['CGCS2000'];
    const e2 = 2*f - f*f;
    y -= 500000;
    
    const e4 = e2*e2, e6 = e4*e2;
    const A0 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
    const A2 = 3/8 * (e2 + e4/4 + 15*e6/128);
    const A4 = 15/256 * (e4 + 3*e6/4);
    const A6 = 35*e6/3072;
    
    // 迭代求底点纬度
    let Bf = x / (a * A0);
    for (let i = 0; i < 10; i++) {
      const FBf = a * (A0*Bf - A2*Math.sin(2*Bf) + A4*Math.sin(4*Bf) - A6*Math.sin(6*Bf)) - x;
      const dF = a * A0 * (1 - e2*Math.sin(Bf)*Math.sin(Bf));
      Bf -= FBf / dF;
    }
    
    const sinBf = Math.sin(Bf), cosBf = Math.cos(Bf);
    const Nf = a / Math.sqrt(1 - e2*sinBf*sinBf);
    const tf = Math.tan(Bf), tf2 = tf*tf;
    const eta2f = e2/(1-e2) * cosBf*cosBf;
    const y2 = y*y, Nf2 = Nf*Nf;
    
    const B = Bf - tf*y2/(2*Nf2)*(1+eta2f) + tf*y2*y2/(24*Nf2*Nf2)*(5+3*tf2+eta2f-9*eta2f*tf2);
    const l = y/(Nf*cosBf) * (1 - y2/(6*Nf2)*(1+2*tf2+eta2f) + y2*y2/(120*Nf2*Nf2)*(5+28*tf2+24*tf2*tf2));
    
    return { B: Survey.radToDeg(B), L: Survey.radToDeg(l) + cm };
  },
  
  // UTM投影
  utm: (B: number, L: number, ellipsoid: string = 'WGS84') => {
    const zone = Math.floor((L + 180) / 6) + 1;
    const cm = zone * 6 - 183;
    const result = Survey.gaussProj(B, L, 6, cm, ellipsoid);
    return {
      x: result.x * 0.9996,
      y: result.y * 0.9996 + (B < 0 ? 10000000 : 0),
      zone,
      cm,
      hemisphere: B >= 0 ? 'N' : 'S'
    };
  },
  
  // 常用坐标系转换参数（近似值）
  transformParams: {
    'WGS84_TO_CGCS2000': { dx: 0, dy: 0, dz: 0, rx: 0, ry: 0, rz: 0, m: 0 }, // 几乎一致
    'WGS84_TO_BJ54': { dx: -12.064, dy: 130.632, dz: 81.99, rx: 1.168, ry: -0.298, rz: 0.301, m: 6.389 },
    'WGS84_TO_XIAN80': { dx: 24, dy: -123, dz: -94, rx: -0.02, ry: 0.353, rz: -0.22, m: -0.9 },
    'CGCS2000_TO_BJ54': { dx: -12, dy: 131, dz: 82, rx: 1.17, ry: -0.3, rz: 0.3, m: 6.4 },
    'CGCS2000_TO_XIAN80': { dx: 24, dy: -123, dz: -94, rx: -0.02, ry: 0.35, rz: -0.22, m: -0.9 }
  } as {[k:string]:{dx:number;dy:number;dz:number;rx:number;ry:number;rz:number;m:number}}
};

// ==================== 主应用 ====================
function App() {
  const [tab, setTab] = useState<'home'|'calc'|'survey'|'settings'|'help'>('home');
  const [display, setDisplay] = useState('0');
  const [expr, setExpr] = useState('');
  const [calcExpr, setCalcExpr] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [historyFilter, setHistoryFilter] = useState<'all'|'calc'|'survey'>('all');
  const [mem, setMem] = useState(0);
  const [hasMem, setHasMem] = useState(false);
  const [angleUnit, setAngleUnit] = useState<'DEG'|'RAD'|'GRAD'>('DEG');
  const [precision, setPrecision] = useState(6);
  const [vibration, setVibration] = useState(true);
  const [theme, setTheme] = useState<string>('dark');
  const inputRef = useRef<HTMLInputElement>(null);
  
  // 测绘计算状态
  const [surveyType, setSurveyType] = useState('forward');
  const [inputs, setInputs] = useState<{[k:string]:string}>({});
  const [result, setResult] = useState('');
  
  // 获取当前主题配色
  const currentTheme = Themes[theme] || Themes['dark'];
  
  useEffect(() => {
    const saved = localStorage.getItem('survey_history');
    if(saved) setHistory(JSON.parse(saved));
    const savedTheme = localStorage.getItem('survey_theme');
    if(savedTheme) setTheme(savedTheme);
  }, []);
  
  // 应用主题色
  useEffect(() => {
    document.documentElement.style.setProperty('--bg-color', currentTheme.bg);
    document.documentElement.style.setProperty('--card-color', currentTheme.card);
    document.documentElement.style.setProperty('--primary-color', currentTheme.primary);
    document.documentElement.style.setProperty('--text-color', currentTheme.text);
    document.documentElement.style.setProperty('--border-color', currentTheme.border);
  }, [currentTheme]);
  
  const vibrate = useCallback(() => {
    if(vibration && navigator.vibrate) navigator.vibrate(15);
  }, [vibration]);
  
  const fmt = (n: number) => n.toFixed(precision);
  
  const saveCalcHistory = (e: string, r: string) => {
    const item: HistoryItem = { expression: e, result: r, time: Date.now(), category: 'calc' };
    const h = [item, ...history].slice(0, 100);
    setHistory(h);
    localStorage.setItem('survey_history', JSON.stringify(h));
  };
  
  const saveSurveyHistory = (e: string, r: string) => {
    const item: HistoryItem = { expression: e, result: r, time: Date.now(), category: 'survey', surveyType, inputs: {...inputs} };
    const h = [item, ...history].slice(0, 100);
    setHistory(h);
    localStorage.setItem('survey_history', JSON.stringify(h));
  };
  
  // 点击历史记录跳转
  const jumpToHistory = (item: HistoryItem) => {
    vibrate();
    if(item.category === 'calc') {
      setTab('calc');
      setDisplay(item.result);
      setExpr(item.expression + ' =');
    } else if(item.category === 'survey' && item.surveyType) {
      setTab('survey');
      setSurveyType(item.surveyType);
      if(item.inputs) setInputs(item.inputs);
      setResult(item.result);
    }
  };
  
  // 切换主题
  const changeTheme = (t: string) => {
    setTheme(t);
    localStorage.setItem('survey_theme', t);
  };
  
  // 计算器函数
  const clear = () => { vibrate(); setDisplay('0'); setExpr(''); setCalcExpr(''); };
  const append = (v: string) => {
    vibrate();
    if(display === '0' && v !== '.') { setDisplay(v); setCalcExpr(v); }
    else if(display === 'Error') { setDisplay(v); setCalcExpr(v); }
    else { setDisplay(display + v); setCalcExpr(calcExpr + v); }
  };
  const back = () => { vibrate(); setDisplay(display.length > 1 ? display.slice(0,-1) : '0'); setCalcExpr(calcExpr.length > 1 ? calcExpr.slice(0,-1) : ''); };
  const toggleSign = () => { vibrate(); if(display !== '0') setDisplay(display.startsWith('-') ? display.slice(1) : '-'+display); };
    
  const calc = () => {
    vibrate();
    try {
      let e = display
        .replace(/×/g,'*').replace(/÷/g,'/').replace(/π/g,`(${Math.PI})`).replace(/\^/g,'**')
        .replace(/√\(/g,'Math.sqrt(').replace(/∛\(/g,'Math.cbrt(')
        .replace(/sin\(/g,`Math.sin(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''})`)
        .replace(/cos\(/g,`Math.cos(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''})`)
        .replace(/tan\(/g,`Math.tan(${angleUnit==='DEG'?'Math.PI/180*':angleUnit==='GRAD'?'Math.PI/200*':''})`)
        .replace(/asin\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.asin()`)
        .replace(/acos\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.acos()`)
        .replace(/atan\(/g,`(${angleUnit==='DEG'?'180/Math.PI*':angleUnit==='GRAD'?'200/Math.PI*':''}Math.atan()`)
        .replace(/ln\(/g,'Math.log(').replace(/log\(/g,'Math.log10(')
        .replace(/abs\(/g,'Math.abs(').replace(/exp\(/g,'Math.exp(');
      const r = eval(e);
      const res = fmt(r);
      saveCalcHistory(display, res);
      setExpr(display + ' =');
      setDisplay(res);
      setCalcExpr('');
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
      const exprStr = `${fn}(${display})`;
      setExpr(exprStr + ' =');
      saveCalcHistory(exprStr, fmt(r));
      setDisplay(fmt(r));
    } catch { setDisplay('Error'); }
  };
  
  const insertFn = (fn: string) => { 
    vibrate(); 
    const newExpr = display==='0'||display==='Error' ? fn+'(' : display+fn+'(';
    setDisplay(newExpr);
    setExpr(newExpr);
  };
  
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
      setExpr(`D→DMS(${display}) =`);
      saveCalcHistory(`D→DMS(${display})`, r);
      setDisplay(r);
    } catch { setDisplay('Error'); }
  };
    
  const dms2deg = () => {
    vibrate();
    try {
      const inp = display.trim();
      let deg = 0;
      if(inp.includes('°')) {
        const parts = inp.replace(/['"\u2032\u2033]/g,' ').replace('°',' ').trim().split(/\s+/);
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
      setExpr(`DMS→D(${display}) =`);
      saveCalcHistory(`DMS→D(${display})`, r);
      setDisplay(r);
    } catch { setDisplay('Error'); }
  };
    
  const deg2rad = () => { vibrate(); try { const r = fmt(parseFloat(display)*Math.PI/180); setExpr(`D→R(${display}) =`); saveCalcHistory(`D→R(${display})`,r); setDisplay(r); } catch { setDisplay('Error'); } };
  const rad2deg = () => { vibrate(); try { const r = fmt(parseFloat(display)*180/Math.PI); setExpr(`R→D(${display}) =`); saveCalcHistory(`R→D(${display})`,r); setDisplay(r); } catch { setDisplay('Error'); } };
  
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
        case 'gauss_proj': {
          const zw = inputs['gpzw'] === '3' ? 3 : 6;
          const ellip = inputs['gpellip'] || 'CGCS2000';
          const g = Survey.gaussProj(getN('gpB'), getN('gpL'), zw as 3|6, getN('gpL0')||undefined, ellip);
          r = '【高斯投影结果】\n\n椭球: ' + (Ellipsoids[ellip]?.name || ellip) + '\n带宽: ' + zw + '°\n\n输入:\nB = ' + fmt(getN('gpB')) + '°\nL = ' + fmt(getN('gpL')) + '°\n\n输出:\nX = ' + fmt(g.x) + ' m\nY = ' + fmt(g.y) + ' m\n带号 = ' + g.zone + '\n中央子午线 = ' + g.cm + '°';
          break;
        }
        case 'utm': {
          const ellip = inputs['utmellip'] || 'WGS84';
          const u = Survey.utm(getN('utmB'), getN('utmL'), ellip);
          r = '【UTM投影结果】\n\n椭球: ' + (Ellipsoids[ellip]?.name || ellip) + '\n\n输入:\nB = ' + fmt(getN('utmB')) + '°\nL = ' + fmt(getN('utmL')) + '°\n\n输出:\nN(X) = ' + fmt(u.x) + ' m\nE(Y) = ' + fmt(u.y) + ' m\n带号 = ' + u.zone + u.hemisphere + '\n中央子午线 = ' + u.cm + '°';
          break;
        }
        case 'blh_xyz': {
          const mode = inputs['blhmode'] || 'blh2xyz';
          const ellip = inputs['blhellip'] || 'CGCS2000';
          if (mode === 'blh2xyz') {
            const xyz = Survey.blhToXyz(getN('blhB'), getN('blhL'), getN('blhH'), ellip);
            r = '【BLH→XYZ转换】\n\n椭球: ' + (Ellipsoids[ellip]?.name || ellip) + '\n\n输入(大地坐标):\nB = ' + fmt(getN('blhB')) + '°\nL = ' + fmt(getN('blhL')) + '°\nH = ' + fmt(getN('blhH')) + ' m\n\n输出(空间直角坐标):\nX = ' + fmt(xyz.X) + ' m\nY = ' + fmt(xyz.Y) + ' m\nZ = ' + fmt(xyz.Z) + ' m';
          } else {
            const blh = Survey.xyzToBlh(getN('xyzX'), getN('xyzY'), getN('xyzZ'), ellip);
            r = '【XYZ→BLH转换】\n\n椭球: ' + (Ellipsoids[ellip]?.name || ellip) + '\n\n输入(空间直角坐标):\nX = ' + fmt(getN('xyzX')) + ' m\nY = ' + fmt(getN('xyzY')) + ' m\nZ = ' + fmt(getN('xyzZ')) + ' m\n\n输出(大地坐标):\nB = ' + fmt(blh.B) + '° (' + Survey.formatDms(blh.B) + ')\nL = ' + fmt(blh.L) + '° (' + Survey.formatDms(blh.L) + ')\nH = ' + fmt(blh.H) + ' m';
          }
          break;
        }
        case 'transform7': {
          const mode = inputs['t7mode'] || 'calc';
          if (mode === 'calc') {
            const src: {X:number;Y:number;Z:number}[] = [], tgt: {X:number;Y:number;Z:number}[] = [];
            for(let i=1; i<=5; i++) {
              const sx=inputs[`t7sX${i}`], sy=inputs[`t7sY${i}`], sz=inputs[`t7sZ${i}`];
              const tx=inputs[`t7tX${i}`], ty=inputs[`t7tY${i}`], tz=inputs[`t7tZ${i}`];
              if(sx&&sy&&sz&&tx&&ty&&tz) {
                src.push({X:parseFloat(sx),Y:parseFloat(sy),Z:parseFloat(sz)});
                tgt.push({X:parseFloat(tx),Y:parseFloat(ty),Z:parseFloat(tz)});
              }
            }
            if(src.length < 3) { r = '至少需要3个公共点'; break; }
            const p = Survey.calc7Param(src, tgt);
            r = '【七参数求解结果】\n\n公共点数: ' + src.length + '\n\n布尔萨参数:\nΔX = ' + fmt(p.dx) + ' m\nΔY = ' + fmt(p.dy) + ' m\nΔZ = ' + fmt(p.dz) + ' m\nεx = ' + fmt(p.rx) + '"\nεy = ' + fmt(p.ry) + '"\nεz = ' + fmt(p.rz) + '"\nm = ' + fmt(p.m) + ' ppm';
          } else {
            const xyz = Survey.transform7Param(getN('t7X'), getN('t7Y'), getN('t7Z'), getN('t7dx'), getN('t7dy'), getN('t7dz'), getN('t7rx'), getN('t7ry'), getN('t7rz'), getN('t7m'));
            r = '【七参数转换结果】\n\n输入:\nX = ' + fmt(getN('t7X')) + ' m\nY = ' + fmt(getN('t7Y')) + ' m\nZ = ' + fmt(getN('t7Z')) + ' m\n\n参数:\nΔX=' + getN('t7dx') + ', ΔY=' + getN('t7dy') + ', ΔZ=' + getN('t7dz') + '\nεx=' + getN('t7rx') + '", εy=' + getN('t7ry') + '", εz=' + getN('t7rz') + '"\nm=' + getN('t7m') + 'ppm\n\n输出:\nX\' = ' + fmt(xyz.X) + ' m\nY\' = ' + fmt(xyz.Y) + ' m\nZ\' = ' + fmt(xyz.Z) + ' m';
          }
          break;
        }
        case 'coord_sys': {
          const srcSys = inputs['csSrc'] || 'WGS84';
          const tgtSys = inputs['csTgt'] || 'CGCS2000';
          const paramKey = srcSys + '_TO_' + tgtSys;
          const params = Survey.transformParams[paramKey];
          
          if (!params) {
            r = '暂不支持' + srcSys + '→' + tgtSys + '转换\n\n支持的转换:\nWGS84→CGCS2000\nWGS84→BJ54\nWGS84→XIAN80\nCGCS2000→BJ54\nCGCS2000→XIAN80';
            break;
          }
          
          const srcXyz = Survey.blhToXyz(getN('csB'), getN('csL'), getN('csH'), srcSys);
          const tgtXyz = Survey.transform7Param(srcXyz.X, srcXyz.Y, srcXyz.Z, params.dx, params.dy, params.dz, params.rx, params.ry, params.rz, params.m);
          const tgtBlh = Survey.xyzToBlh(tgtXyz.X, tgtXyz.Y, tgtXyz.Z, tgtSys);
          
          r = '【坐标系转换结果】\n\n' + srcSys + ' → ' + tgtSys + '\n\n源坐标:' + '\nB = ' + fmt(getN('csB')) + '°\nL = ' + fmt(getN('csL')) + '°\nH = ' + fmt(getN('csH')) + ' m\n\n转换参数:\nΔX=' + params.dx + 'm, ΔY=' + params.dy + 'm, ΔZ=' + params.dz + 'm\n\n目标坐标:\nB = ' + fmt(tgtBlh.B) + '° (' + Survey.formatDms(tgtBlh.B) + ')\nL = ' + fmt(tgtBlh.L) + '° (' + Survey.formatDms(tgtBlh.L) + ')\nH = ' + fmt(tgtBlh.H) + ' m';
          break;
        }
        default: r = '请选择计算类型';
      }
      setResult(r);
      saveSurveyHistory(surveyType, r.split('\n')[0]);
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
    { id: 'gauss_proj', name: '高斯投影(3°/6°)', icon: '📜' },
    { id: 'utm', name: 'UTM投影', icon: '🌐' },
    { id: 'blh_xyz', name: 'BLH↔XYZ', icon: '📡' },
    { id: 'transform4', name: '四参数转换', icon: '🔄' },
    { id: 'transform7', name: '七参数转换', icon: '🔀' },
    { id: 'coord_sys', name: '坐标系转换', icon: '🔁' },
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
      case 'gauss_proj':
        return <><div className="select-row"><label>椭球</label><select value={inputs['gpellip']||'CGCS2000'} onChange={e=>inp('gpellip',e.target.value)}><option value="CGCS2000">CGCS2000</option><option value="WGS84">WGS84</option><option value="BJ54">北京54</option><option value="XIAN80">西安80</option></select></div><div className="select-row"><label>带宽</label><select value={inputs['gpzw']||'6'} onChange={e=>inp('gpzw',e.target.value)}><option value="6">6°带</option><option value="3">3°带</option></select></div><InputField label="纬度B(°)" k="gpB" placeholder="如 30.5"/><InputField label="经度L(°)" k="gpL" placeholder="如 114.3"/><InputField label="中央子午线(°)" k="gpL0" placeholder="自动计算"/></>;
      case 'utm':
        return <><div className="select-row"><label>椭球</label><select value={inputs['utmellip']||'WGS84'} onChange={e=>inp('utmellip',e.target.value)}><option value="WGS84">WGS84</option><option value="CGCS2000">CGCS2000</option></select></div><InputField label="纬度B(°)" k="utmB" placeholder="如 30.5"/><InputField label="经度L(°)" k="utmL" placeholder="如 114.3"/></>;
      case 'blh_xyz':
        return <><div className="select-row"><label>椭球</label><select value={inputs['blhellip']||'CGCS2000'} onChange={e=>inp('blhellip',e.target.value)}><option value="CGCS2000">CGCS2000</option><option value="WGS84">WGS84</option><option value="BJ54">北京54</option><option value="XIAN80">西安80</option></select></div><div className="select-row"><label>转换方向</label><select value={inputs['blhmode']||'blh2xyz'} onChange={e=>inp('blhmode',e.target.value)}><option value="blh2xyz">BLH→XYZ</option><option value="xyz2blh">XYZ→BLH</option></select></div>{(inputs['blhmode']||'blh2xyz')==='blh2xyz'?<><InputField label="纬度B(°)" k="blhB"/><InputField label="经度L(°)" k="blhL"/><InputField label="大地高H(m)" k="blhH"/></>:<><InputField label="X(m)" k="xyzX"/><InputField label="Y(m)" k="xyzY"/><InputField label="Z(m)" k="xyzZ"/></>}</>;
      case 'transform7':
        return <><div className="select-row"><label>模式</label><select value={inputs['t7mode']||'calc'} onChange={e=>inp('t7mode',e.target.value)}><option value="calc">参数求解</option><option value="apply">参数转换</option></select></div>{(inputs['t7mode']||'calc')==='calc'?<><div className="transform-header">公共点坐标（至少3个）- 空间直角坐标</div><div className="table-header"><span>点</span><span>源X</span><span>源Y</span><span>源Z</span></div>{[1,2,3,4,5].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`t7sX${i}`]||''} onChange={e=>inp(`t7sX${i}`,e.target.value)} placeholder="X"/><input type="text" inputMode="decimal" value={inputs[`t7sY${i}`]||''} onChange={e=>inp(`t7sY${i}`,e.target.value)} placeholder="Y"/><input type="text" inputMode="decimal" value={inputs[`t7sZ${i}`]||''} onChange={e=>inp(`t7sZ${i}`,e.target.value)} placeholder="Z"/></div>)}<div className="table-header"><span>点</span><span>目X</span><span>目Y</span><span>目Z</span></div>{[1,2,3,4,5].map(i=><div key={i} className="table-row"><span>{i}</span><input type="text" inputMode="decimal" value={inputs[`t7tX${i}`]||''} onChange={e=>inp(`t7tX${i}`,e.target.value)} placeholder="X'"/><input type="text" inputMode="decimal" value={inputs[`t7tY${i}`]||''} onChange={e=>inp(`t7tY${i}`,e.target.value)} placeholder="Y'"/><input type="text" inputMode="decimal" value={inputs[`t7tZ${i}`]||''} onChange={e=>inp(`t7tZ${i}`,e.target.value)} placeholder="Z'"/></div>)}</>:<><div className="transform-header">布尔萨七参数</div><InputField label="ΔX(m)" k="t7dx"/><InputField label="ΔY(m)" k="t7dy"/><InputField label="ΔZ(m)" k="t7dz"/><InputField label="εx(角秒)" k="t7rx"/><InputField label="εy(角秒)" k="t7ry"/><InputField label="εz(角秒)" k="t7rz"/><InputField label="m(ppm)" k="t7m"/><div className="transform-header">待转换点</div><InputField label="X(m)" k="t7X"/><InputField label="Y(m)" k="t7Y"/><InputField label="Z(m)" k="t7Z"/></>}</>;
      case 'coord_sys':
        return <><div className="select-row"><label>源坐标系</label><select value={inputs['csSrc']||'WGS84'} onChange={e=>inp('csSrc',e.target.value)}><option value="WGS84">WGS84</option><option value="CGCS2000">CGCS2000</option></select></div><div className="select-row"><label>目标坐标系</label><select value={inputs['csTgt']||'CGCS2000'} onChange={e=>inp('csTgt',e.target.value)}><option value="CGCS2000">CGCS2000</option><option value="BJ54">北京54</option><option value="XIAN80">西安80</option></select></div><InputField label="纬度B(°)" k="csB"/><InputField label="经度L(°)" k="csL"/><InputField label="大地高H(m)" k="csH"/></>;
      default: return null;
    }
  };

  return (
    <div className={`app ${theme}`} style={{background: currentTheme.bg, color: currentTheme.text}}>
      <div className="main-content">
        {/* 首页 */}
        {tab === 'home' && (
          <div className="home-page">
            <h1>测绘计算器Pro</h1>
            <p className="subtitle">专业测绘计算 · 向导式操作</p>
            <div className="quick-grid">
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('calc');}}>
                <span className="icon">📱</span>
                <span>科学计算</span>
              </div>
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('survey');setSurveyType('forward');}}>
                <span className="icon">📍</span>
                <span>坐标正反算</span>
              </div>
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('survey');setSurveyType('closed_traverse');}}>
                <span className="icon">🔄</span>
                <span>导线计算</span>
              </div>
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('survey');setSurveyType('level_closed');}}>
                <span className="icon">📊</span>
                <span>水准平差</span>
              </div>
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('survey');setSurveyType('gauss_proj');}}>
                <span className="icon">🌍</span>
                <span>高斯投影</span>
              </div>
              <div className="quick-card" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>{setTab('survey');setSurveyType('coord_sys');}}>
                <span className="icon">🔀</span>
                <span>坐标转换</span>
              </div>
            </div>
            <div className="history-section">
              <div className="history-header">
                <h3>历史记录</h3>
                <div className="history-filter">
                  <button className={historyFilter==='all'?'active':''} onClick={()=>setHistoryFilter('all')}>全部</button>
                  <button className={historyFilter==='calc'?'active':''} onClick={()=>setHistoryFilter('calc')}>计算器</button>
                  <button className={historyFilter==='survey'?'active':''} onClick={()=>setHistoryFilter('survey')}>测绘</button>
                </div>
              </div>
              {history.filter(h => historyFilter === 'all' || h.category === historyFilter).length === 0 ? <p className="empty">暂无记录</p> : (
                <div className="history-list">
                  {history.filter(h => historyFilter === 'all' || h.category === historyFilter).slice(0,15).map((h,i) => (
                    <div key={i} className="history-item" onClick={()=>jumpToHistory(h)} style={{background: currentTheme.card, borderColor: currentTheme.border}}>
                      <div className="history-item-left">
                        <span className="history-category" style={{background: h.category==='calc' ? '#1f6feb' : currentTheme.primary}}>{h.category==='calc'?'计算':'测绘'}</span>
                        <span className="expr">{h.expression}</span>
                      </div>
                      <span className="res" style={{color: currentTheme.primary}}>{h.result.length > 20 ? h.result.substring(0,20)+'...' : h.result}</span>
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
            <div className="display-area" style={{background: currentTheme.card, borderColor: currentTheme.border}}>
              <div className="expression" style={{color: currentTheme.text, opacity: 0.6}}>{expr}</div>
              <div className="display" style={{color: currentTheme.text}}>{display}</div>
              <div className="status-bar">
                <span className={hasMem ? 'active' : ''} style={{color: hasMem ? currentTheme.primary : undefined}}>M</span>
                <span style={{color: currentTheme.primary}}>{angleUnit}</span>
              </div>
            </div>
            <div className="sci-panel">
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('sin')}>sin</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('cos')}>cos</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('tan')}>tan</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('sinh')}>sinh</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('cosh')}>cosh</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('asin')}>sin⁻¹</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('acos')}>cos⁻¹</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('atan')}>tan⁻¹</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('ln')}>ln</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('log')}>log</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('√')}>√</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('∛')}>∛</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('x²')}>x²</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('x³')}>x³</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append('^')}>^</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('n!')}>n!</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('1/x')}>1/x</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('abs')}>|x|</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('eˣ')}>eˣ</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('10ˣ')}>10ˣ</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={deg2dms}>D→DMS</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={dms2deg}>DMS→D</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={deg2rad}>D→R</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={rad2deg}>R→D</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append('°')}>°</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={mc}>MC</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={mr}>MR</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={mAdd}>M+</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={mSub}>M-</button>
                <button style={{background: currentTheme.primary, color: '#fff', borderColor: currentTheme.border}} onClick={()=>setAngleUnit(angleUnit==='DEG'?'RAD':angleUnit==='RAD'?'GRAD':'DEG')}>{angleUnit}</button>
              </div>
              <div className="sci-row">
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append('π')}>π</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append('e')}>e</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append('(')}>{'('}</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>append(')')}>{')'}</button>
                <button style={{background: currentTheme.card, color: currentTheme.text, borderColor: currentTheme.border}} onClick={()=>applyFn('%')}>%</button>
              </div>
            </div>
            <div className="num-panel">
              <div className="num-row">
                <button className="func" style={{background: currentTheme.card, color: '#f78166'}} onClick={clear}>C</button>
                <button className="func" style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>setDisplay('0')}>CE</button>
                <button className="func" style={{background: currentTheme.card, color: currentTheme.text}} onClick={back}>⌫</button>
                <button className="op" style={{background: currentTheme.primary}} onClick={()=>append('÷')}>÷</button>
              </div>
              <div className="num-row">
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('7')}>7</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('8')}>8</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('9')}>9</button>
                <button className="op" style={{background: currentTheme.primary}} onClick={()=>append('×')}>×</button>
              </div>
              <div className="num-row">
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('4')}>4</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('5')}>5</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('6')}>6</button>
                <button className="op" style={{background: currentTheme.primary}} onClick={()=>append('-')}>-</button>
              </div>
              <div className="num-row">
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('1')}>1</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('2')}>2</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('3')}>3</button>
                <button className="op" style={{background: currentTheme.primary}} onClick={()=>append('+')}>+</button>
              </div>
              <div className="num-row">
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={toggleSign}>±</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('0')}>0</button>
                <button style={{background: currentTheme.card, color: currentTheme.text}} onClick={()=>append('.')}>.</button>
                <button className="eq" style={{background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={calc}>=</button>
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
              <button className={`toggle ${vibration?'on':''}`} style={{background: vibration ? currentTheme.primary : currentTheme.card}} onClick={()=>setVibration(!vibration)}>{vibration?'开':'关'}</button>
            </div>
            <div className="setting-item" style={{borderColor: currentTheme.border}}>
              <span>主题色彩</span>
              <div className="theme-selector">
                {Object.entries(Themes).map(([key, t]) => (
                  <button 
                    key={key} 
                    className={`theme-btn ${theme===key?'active':''}`}
                    style={{background: t.primary, border: theme===key ? '3px solid #fff' : 'none'}}
                    onClick={()=>changeTheme(key)}
                    title={t.name}
                  >{t.name.charAt(0)}</button>
                ))}
              </div>
            </div>
            <div className="setting-item" style={{borderColor: currentTheme.border}}>
              <span>清除历史</span>
              <button className="danger" onClick={()=>{setHistory([]);localStorage.removeItem('survey_history');}}>清除</button>
            </div>
            <div className="setting-item" style={{borderColor: currentTheme.border}}>
              <span>帮助与反馈</span>
              <button className="toggle on" style={{background: currentTheme.primary}} onClick={()=>setTab('help')}>查看</button>
            </div>
            <div className="about">
              <p>测绘计算器Pro v3.0</p>
              <p>专业测绘计算解决方案</p>
            </div>
          </div>
        )}
        
        {/* 帮助页面 */}
        {tab === 'help' && (
          <div className="help-page" style={{padding: 20}}>
            <h2>帮助与知识库</h2>
            
            <div className="help-section" style={{marginTop: 20}}>
              <h3>📚 坐标转换知识</h3>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <p><strong>七参数转换</strong>：适用于不同坐标系间的精确转换，包含3个平移(ΔX/ΔY/ΔZ)、3个旋转(εx/εy/εz)和1个尺度因子(m).</p>
              </div>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <p><strong>四参数转换</strong>：适用于小范围平面坐标转换，包含2个平移、1个旋转和1个尺度变化.</p>
              </div>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <p><strong>高斯投影</strong>：将大地坐标(BLH)投影到平面坐标，支持3°带和6°带.</p>
              </div>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <p><strong>UTM投影</strong>：通用横轴墨卡托投影，将地球划分60个投影带.</p>
              </div>
            </div>
            
            <div className="help-section" style={{marginTop: 24}}>
              <h3>🌐 支持的坐标系</h3>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <ul style={{paddingLeft: 20, lineHeight: 2}}>
                  <li><strong>CGCS2000</strong> - 2000国家大地坐标系</li>
                  <li><strong>WGS84</strong> - GPS全球定位系统</li>
                  <li><strong>北京54</strong> - 1954年北京坐标系</li>
                  <li><strong>西安80</strong> - 1980年西安坐标系</li>
                </ul>
              </div>
            </div>
            
            <div className="help-section" style={{marginTop: 24}}>
              <h3>🔗 相关链接</h3>
              <div className="help-card" style={{background: currentTheme.card, padding: 16, borderRadius: 12, marginTop: 12}}>
                <p onClick={()=>window.open('https://github.com/mzy222603/survey-calculator-app')} style={{color: currentTheme.primary, cursor: 'pointer'}}>📦 GitHub仓库 - 查看源码和更新</p>
              </div>
            </div>
            
            <button className="calc-btn" style={{marginTop: 24, background: `linear-gradient(135deg, ${currentTheme.primary}, #1f6feb)`}} onClick={()=>setTab('settings')}>返回设置</button>
          </div>
        )}
      </div>

      {/* 底部导航 */}
      <nav className="bottom-nav" style={{background: currentTheme.card, borderColor: currentTheme.border}}>
        <button className={tab==='home'?'active':''} style={{color: tab==='home' ? currentTheme.primary : undefined}} onClick={()=>setTab('home')}><span>🏠</span>首页</button>
        <button className={tab==='calc'?'active':''} style={{color: tab==='calc' ? currentTheme.primary : undefined}} onClick={()=>setTab('calc')}><span>📱</span>计算器</button>
        <button className={tab==='survey'?'active':''} style={{color: tab==='survey' ? currentTheme.primary : undefined}} onClick={()=>setTab('survey')}><span>📐</span>测绘</button>
        <button className={tab==='settings'||tab==='help'?'active':''} style={{color: (tab==='settings'||tab==='help') ? currentTheme.primary : undefined}} onClick={()=>setTab('settings')}><span>⚙️</span>设置</button>
      </nav>
    </div>
  );
}

export default App;
