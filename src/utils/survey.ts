/**
 * 测绘计算引擎
 * 专业测量计算功能
 */

export interface Point {
  x: number;  // 纵坐标 N
  y: number;  // 横坐标 E
  z?: number; // 高程
  name?: string;
}

export interface TraverseStation {
  angle: number;      // 观测角（度）
  distance: number;   // 边长（米）
  heightDiff?: number; // 高差
}

export interface TraverseResult {
  points: Point[];
  angleClosure: number;
  coordinateClosure: { fx: number; fy: number; f: number };
  relativeClosure: number;
  adjustedAngles: number[];
  adjustedDistances: { dx: number; dy: number }[];
}

export interface GaussResult {
  x: number;  // 纵坐标
  y: number;  // 横坐标
  zone: number; // 带号
  centralMeridian: number; // 中央子午线
}

export interface CurveElements {
  radius: number;
  deflectionAngle: number;
  tangentLength: number;
  curveLength: number;
  externalDistance: number;
  chord: number;
}

// 椭球参数
const ELLIPSOID = {
  CGCS2000: { a: 6378137.0, f: 1/298.257222101 },
  WGS84: { a: 6378137.0, f: 1/298.257223563 },
  BJ54: { a: 6378245.0, f: 1/298.3 },
  XIAN80: { a: 6378140.0, f: 1/298.257 },
};

class SurveyCalculator {
  private ellipsoid = ELLIPSOID.CGCS2000;

  setEllipsoid(name: keyof typeof ELLIPSOID): void {
    this.ellipsoid = ELLIPSOID[name];
  }

  // ==================== 角度转换 ====================
  
  degToRad(deg: number): number { return deg * Math.PI / 180; }
  radToDeg(rad: number): number { return rad * 180 / Math.PI; }
  
  dmsToD(d: number, m: number, s: number): number {
    const sign = d >= 0 ? 1 : -1;
    return sign * (Math.abs(d) + m / 60 + s / 3600);
  }
  
  dToDms(deg: number): { d: number; m: number; s: number } {
    const sign = deg >= 0 ? 1 : -1;
    deg = Math.abs(deg);
    const d = Math.floor(deg);
    const mFull = (deg - d) * 60;
    const m = Math.floor(mFull);
    const s = (mFull - m) * 60;
    return { d: sign * d, m, s };
  }
  
  formatDms(deg: number): string {
    const { d, m, s } = this.dToDms(deg);
    return `${d}°${m}'${s.toFixed(2)}"`;
  }
  
  normalizeAzimuth(azimuth: number): number {
    while (azimuth < 0) azimuth += 360;
    while (azimuth >= 360) azimuth -= 360;
    return azimuth;
  }

  // ==================== 坐标计算 ====================
  
  // 坐标正算
  forwardCalc(point: Point, azimuth: number, distance: number): Point {
    const azRad = this.degToRad(azimuth);
    return {
      x: point.x + distance * Math.cos(azRad),
      y: point.y + distance * Math.sin(azRad),
      z: point.z
    };
  }
  
  // 坐标反算
  inverseCalc(p1: Point, p2: Point): { azimuth: number; distance: number } {
    const dx = p2.x - p1.x;
    const dy = p2.y - p1.y;
    const distance = Math.sqrt(dx * dx + dy * dy);
    let azimuth = this.radToDeg(Math.atan2(dy, dx));
    azimuth = this.normalizeAzimuth(azimuth);
    return { azimuth, distance };
  }
  
  // 前方交会
  forwardIntersection(pa: Point, pb: Point, angleA: number, angleB: number): Point {
    const { azimuth: azAB } = this.inverseCalc(pa, pb);
    const azAP = this.normalizeAzimuth(azAB - angleA);
    const azBP = this.normalizeAzimuth(azAB + 180 + angleB);
    
    const alphaRad = this.degToRad(azAP);
    const betaRad = this.degToRad(azBP);
    
    const sinA = Math.sin(alphaRad);
    const cosA = Math.cos(alphaRad);
    const sinB = Math.sin(betaRad);
    const cosB = Math.cos(betaRad);
    
    const denom = sinA * cosB - cosA * sinB;
    if (Math.abs(denom) < 1e-10) throw new Error('两条方向线平行');
    
    const t = ((pb.y - pa.y) * cosB - (pb.x - pa.x) * sinB) / denom;
    
    return {
      x: pa.x + t * cosA,
      y: pa.y + t * sinA
    };
  }
  
  // 后方交会（三点法）
  resection(pa: Point, pb: Point, pc: Point, alphaDeg: number, betaDeg: number): Point {
    // 使用科特公式
    const alpha = this.degToRad(alphaDeg);
    const beta = this.degToRad(betaDeg);
    
    const cotA = 1 / Math.tan(alpha);
    const cotB = 1 / Math.tan(beta);
    
    const k1 = (pb.x - pa.x) * cotA - (pb.y - pa.y);
    const k2 = (pb.y - pa.y) * cotA + (pb.x - pa.x);
    const k3 = (pc.x - pb.x) * cotB - (pc.y - pb.y);
    const k4 = (pc.y - pb.y) * cotB + (pc.x - pb.x);
    
    const denom = k2 - k4;
    if (Math.abs(denom) < 1e-10) throw new Error('无法求解');
    
    const y = (k1 - k3) / denom;
    const x = pa.x + (pa.y - y) * cotA;
    
    return { x, y };
  }
  
  // 侧方交会
  sideIntersection(pa: Point, pb: Point, angleA: number, distAP: number): Point {
    const { azimuth: azAB } = this.inverseCalc(pa, pb);
    const azAP = this.normalizeAzimuth(azAB + angleA);
    return this.forwardCalc(pa, azAP, distAP);
  }

  // ==================== 导线计算 ====================
  
  // 闭合导线计算
  closedTraverse(startPoint: Point, startAzimuth: number, stations: TraverseStation[]): TraverseResult {
    const n = stations.length;
    if (n < 3) throw new Error('闭合导线至少需要3个测站');
    
    // 理论角度和
    const theoreticalSum = (n - 2) * 180;
    const measuredSum = stations.reduce((sum, s) => sum + s.angle, 0);
    const angleClosure = measuredSum - theoreticalSum;
    
    // 角度改正
    const angleCorrection = -angleClosure / n;
    const adjustedAngles = stations.map(s => s.angle + angleCorrection);
    
    // 推算方位角
    const azimuths: number[] = [startAzimuth];
    for (let i = 0; i < n; i++) {
      const newAz = this.normalizeAzimuth(azimuths[i] + adjustedAngles[i] - 180);
      azimuths.push(newAz);
    }
    
    // 计算坐标增量
    const dxList: number[] = [];
    const dyList: number[] = [];
    for (let i = 0; i < n; i++) {
      const azRad = this.degToRad(azimuths[i + 1]);
      dxList.push(stations[i].distance * Math.cos(azRad));
      dyList.push(stations[i].distance * Math.sin(azRad));
    }
    
    // 坐标闭合差
    const fx = dxList.reduce((a, b) => a + b, 0);
    const fy = dyList.reduce((a, b) => a + b, 0);
    const f = Math.sqrt(fx * fx + fy * fy);
    
    const totalLength = stations.reduce((sum, s) => sum + s.distance, 0);
    const relativeClosure = f / totalLength;
    
    // 坐标改正
    const adjustedDistances: { dx: number; dy: number }[] = [];
    for (let i = 0; i < n; i++) {
      const vx = -fx * stations[i].distance / totalLength;
      const vy = -fy * stations[i].distance / totalLength;
      adjustedDistances.push({
        dx: dxList[i] + vx,
        dy: dyList[i] + vy
      });
    }
    
    // 计算各点坐标
    const points: Point[] = [startPoint];
    for (let i = 0; i < n; i++) {
      const prev = points[points.length - 1];
      points.push({
        x: prev.x + adjustedDistances[i].dx,
        y: prev.y + adjustedDistances[i].dy,
        name: `P${i + 1}`
      });
    }
    
    return {
      points,
      angleClosure,
      coordinateClosure: { fx, fy, f },
      relativeClosure,
      adjustedAngles,
      adjustedDistances
    };
  }
  
  // 附合导线计算
  attachedTraverse(
    startPoint: Point, endPoint: Point,
    startAzimuth: number, endAzimuth: number,
    stations: TraverseStation[]
  ): TraverseResult {
    const n = stations.length;
    
    // 推算终点方位角
    let calcAzimuth = startAzimuth;
    for (const s of stations) {
      calcAzimuth = this.normalizeAzimuth(calcAzimuth + s.angle - 180);
    }
    
    // 角度闭合差
    const angleClosure = calcAzimuth - endAzimuth;
    const angleCorrection = -angleClosure / n;
    
    const adjustedAngles = stations.map(s => s.angle + angleCorrection);
    
    // 推算改正后的方位角
    const azimuths: number[] = [startAzimuth];
    for (let i = 0; i < n; i++) {
      const newAz = this.normalizeAzimuth(azimuths[i] + adjustedAngles[i] - 180);
      azimuths.push(newAz);
    }
    
    // 坐标增量
    const dxList: number[] = [];
    const dyList: number[] = [];
    for (let i = 0; i < n; i++) {
      const azRad = this.degToRad(azimuths[i + 1]);
      dxList.push(stations[i].distance * Math.cos(azRad));
      dyList.push(stations[i].distance * Math.sin(azRad));
    }
    
    // 坐标闭合差
    const sumDx = dxList.reduce((a, b) => a + b, 0);
    const sumDy = dyList.reduce((a, b) => a + b, 0);
    const actualDx = endPoint.x - startPoint.x;
    const actualDy = endPoint.y - startPoint.y;
    const fx = sumDx - actualDx;
    const fy = sumDy - actualDy;
    const f = Math.sqrt(fx * fx + fy * fy);
    
    const totalLength = stations.reduce((sum, s) => sum + s.distance, 0);
    const relativeClosure = f / totalLength;
    
    // 坐标改正
    const adjustedDistances: { dx: number; dy: number }[] = [];
    let cumDist = 0;
    for (let i = 0; i < n; i++) {
      cumDist += stations[i].distance;
      const vx = -fx * cumDist / totalLength;
      const vy = -fy * cumDist / totalLength;
      adjustedDistances.push({
        dx: dxList[i] + (i === 0 ? vx : vx - (-fx * (cumDist - stations[i].distance) / totalLength)),
        dy: dyList[i] + (i === 0 ? vy : vy - (-fy * (cumDist - stations[i].distance) / totalLength))
      });
    }
    
    // 计算各点坐标
    const points: Point[] = [startPoint];
    let curX = startPoint.x;
    let curY = startPoint.y;
    for (let i = 0; i < n; i++) {
      curX += adjustedDistances[i].dx;
      curY += adjustedDistances[i].dy;
      points.push({ x: curX, y: curY, name: `P${i + 1}` });
    }
    
    return {
      points,
      angleClosure,
      coordinateClosure: { fx, fy, f },
      relativeClosure,
      adjustedAngles,
      adjustedDistances
    };
  }

  // ==================== 面积计算 ====================
  
  polygonArea(points: Point[]): number {
    const n = points.length;
    if (n < 3) throw new Error('至少需要3个顶点');
    
    let area = 0;
    for (let i = 0; i < n; i++) {
      const j = (i + 1) % n;
      area += points[i].x * points[j].y;
      area -= points[j].x * points[i].y;
    }
    return Math.abs(area) / 2;
  }
  
  triangleArea(p1: Point, p2: Point, p3: Point): number {
    return this.polygonArea([p1, p2, p3]);
  }
  
  triangleAreaBySides(a: number, b: number, c: number): number {
    const s = (a + b + c) / 2;
    const sq = s * (s - a) * (s - b) * (s - c);
    if (sq < 0) throw new Error('无法构成三角形');
    return Math.sqrt(sq);
  }

  // ==================== 高斯投影 ====================
  
  // 高斯正算（经纬度→平面坐标）
  gaussForward(latDeg: number, lonDeg: number, centralMeridian?: number): GaussResult {
    const { a, f } = this.ellipsoid;
    const e2 = 2 * f - f * f; // 第一偏心率平方
    const e12 = e2 / (1 - e2);  // 第二偏心率平方
    
    // 自动计算中央子午线（6度带）
    const zone = centralMeridian ? Math.round((centralMeridian + 3) / 6) : Math.floor(lonDeg / 6) + 1;
    const L0 = centralMeridian || zone * 6 - 3;
    
    const B = this.degToRad(latDeg);
    const l = this.degToRad(lonDeg - L0);
    
    // 子午线弧长
    const e4 = e2 * e2;
    const e6 = e4 * e2;
    const A0 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
    const A2 = 3/8 * (e2 + e4/4 + 15*e6/128);
    const A4 = 15/256 * (e4 + 3*e6/4);
    const A6 = 35*e6/3072;
    
    const X = a * (A0*B - A2*Math.sin(2*B) + A4*Math.sin(4*B) - A6*Math.sin(6*B));
    
    // 计算坐标
    const N = a / Math.sqrt(1 - e2 * Math.sin(B) * Math.sin(B));
    const t = Math.tan(B);
    const t2 = t * t;
    const eta2 = e12 * Math.cos(B) * Math.cos(B);
    const cosB = Math.cos(B);
    const l2 = l * l;
    
    const x = X + N * t * cosB * cosB * l2 / 2 * (1 + (5 - t2 + 9*eta2 + 4*eta2*eta2) * cosB*cosB*l2/12 + (61 - 58*t2 + t2*t2) * Math.pow(cosB*l, 4)/360);
    
    let y = N * cosB * l * (1 + (1 - t2 + eta2) * cosB*cosB*l2/6 + (5 - 18*t2 + t2*t2 + 14*eta2 - 58*eta2*t2) * Math.pow(cosB*l, 4)/120);
    
    // 加500km偏移
    y += 500000;
    
    return { x, y, zone, centralMeridian: L0 };
  }
  
  // 高斯反算（平面坐标→经纬度）
  gaussInverse(x: number, y: number, centralMeridian: number): { lat: number; lon: number } {
    const { a, f } = this.ellipsoid;
    const e2 = 2 * f - f * f;
    const e12 = e2 / (1 - e2);
    
    // 去掉500km偏移
    y = y - 500000;
    
    // 迭代求底点纬度
    const e4 = e2 * e2;
    const e6 = e4 * e2;
    const A0 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
    
    let Bf = x / (a * A0);
    for (let i = 0; i < 10; i++) {
      const A2 = 3/8 * (e2 + e4/4 + 15*e6/128);
      const A4 = 15/256 * (e4 + 3*e6/4);
      const A6 = 35*e6/3072;
      const FBf = a * (A0*Bf - A2*Math.sin(2*Bf) + A4*Math.sin(4*Bf) - A6*Math.sin(6*Bf)) - x;
      const dFBf = a * A0 * (1 - e2*Math.sin(Bf)*Math.sin(Bf));
      const newBf = Bf - FBf / dFBf;
      if (Math.abs(newBf - Bf) < 1e-12) break;
      Bf = newBf;
    }
    
    const Nf = a / Math.sqrt(1 - e2 * Math.sin(Bf) * Math.sin(Bf));
    const tf = Math.tan(Bf);
    const tf2 = tf * tf;
    const eta2f = e12 * Math.cos(Bf) * Math.cos(Bf);
    const cosBf = Math.cos(Bf);
    
    const y2 = y * y;
    const Nf2 = Nf * Nf;
    
    const B = Bf - tf * y2 / (2 * Nf2) * (1 + eta2f) + tf * y2*y2 / (24 * Nf2*Nf2) * (5 + 3*tf2 + eta2f - 9*eta2f*tf2);
    const l = y / (Nf * cosBf) * (1 - y2 / (6*Nf2) * (1 + 2*tf2 + eta2f) + y2*y2 / (120*Nf2*Nf2) * (5 + 28*tf2 + 24*tf2*tf2));
    
    return {
      lat: this.radToDeg(B),
      lon: this.radToDeg(l) + centralMeridian
    };
  }

  // ==================== 坐标转换 ====================
  
  // 四参数转换
  transform4Param(points: Point[], dx: number, dy: number, scale: number, rotation: number): Point[] {
    const cosR = Math.cos(rotation);
    const sinR = Math.sin(rotation);
    
    return points.map(p => ({
      x: dx + scale * (p.x * cosR - p.y * sinR),
      y: dy + scale * (p.x * sinR + p.y * cosR),
      z: p.z,
      name: p.name
    }));
  }
  
  // 计算四参数
  calc4Param(sourcePoints: Point[], targetPoints: Point[]): { dx: number; dy: number; scale: number; rotation: number } {
    const n = sourcePoints.length;
    if (n < 2) throw new Error('至少需要2个公共点');
    
    let sumXs = 0, sumYs = 0, sumXt = 0, sumYt = 0;
    let sumXs2 = 0, sumYs2 = 0;
    let sumXsXt = 0, sumYsYt = 0, sumXsYt = 0, sumYsXt = 0;
    
    for (let i = 0; i < n; i++) {
      sumXs += sourcePoints[i].x;
      sumYs += sourcePoints[i].y;
      sumXt += targetPoints[i].x;
      sumYt += targetPoints[i].y;
      sumXs2 += sourcePoints[i].x * sourcePoints[i].x;
      sumYs2 += sourcePoints[i].y * sourcePoints[i].y;
      sumXsXt += sourcePoints[i].x * targetPoints[i].x;
      sumYsYt += sourcePoints[i].y * targetPoints[i].y;
      sumXsYt += sourcePoints[i].x * targetPoints[i].y;
      sumYsXt += sourcePoints[i].y * targetPoints[i].x;
    }
    
    const A = sumXs2 + sumYs2;
    const a = (sumXsXt + sumYsYt) / A;
    const b = (sumYsXt - sumXsYt) / A;
    
    const dx = (sumXt - a * sumXs + b * sumYs) / n;
    const dy = (sumYt - a * sumYs - b * sumXs) / n;
    
    const scale = Math.sqrt(a * a + b * b);
    const rotation = Math.atan2(b, a);
    
    return { dx, dy, scale, rotation };
  }

  // ==================== 曲线计算 ====================
  
  // 圆曲线要素
  circularCurve(radius: number, deflectionAngle: number): CurveElements {
    const alpha = this.degToRad(Math.abs(deflectionAngle));
    
    const T = radius * Math.tan(alpha / 2); // 切线长
    const L = radius * alpha; // 曲线长
    const E = radius * (1 / Math.cos(alpha / 2) - 1); // 外矢距
    const C = 2 * radius * Math.sin(alpha / 2); // 弦长
    
    return {
      radius,
      deflectionAngle,
      tangentLength: T,
      curveLength: L,
      externalDistance: E,
      chord: C
    };
  }
  
  // 缓和曲线要素
  transitionCurve(radius: number, ls: number): { transitionAngle: number; shift: number; tangentIncrement: number; endX: number; endY: number } {
    const beta = ls / (2 * radius);
    const betaDeg = this.radToDeg(beta);
    
    const p = ls * ls / (24 * radius) - Math.pow(ls, 4) / (2688 * Math.pow(radius, 3));
    const m = ls / 2 - Math.pow(ls, 3) / (240 * radius * radius);
    
    const x = ls - Math.pow(ls, 5) / (40 * radius * radius * ls * ls);
    const y = ls * ls / (6 * radius) - Math.pow(ls, 4) / (336 * Math.pow(radius, 3));
    
    return {
      transitionAngle: betaDeg,
      shift: p,
      tangentIncrement: m,
      endX: x,
      endY: y
    };
  }

  // ==================== 水准测量 ====================
  
  // 闭合水准路线平差
  levelClosedAdjustment(knownHeight: number, heightDiffs: number[], distances: number[]): { heights: number[]; closure: number } {
    const n = heightDiffs.length;
    const closure = heightDiffs.reduce((a, b) => a + b, 0);
    const totalDist = distances.reduce((a, b) => a + b, 0);
    
    const heights: number[] = [knownHeight];
    for (let i = 0; i < n; i++) {
      const correction = -closure * distances[i] / totalDist;
      const correctedDiff = heightDiffs[i] + correction;
      heights.push(heights[heights.length - 1] + correctedDiff);
    }
    
    return { heights, closure };
  }
  
  // 附合水准路线平差
  levelAttachedAdjustment(startHeight: number, endHeight: number, heightDiffs: number[], distances: number[]): { heights: number[]; closure: number } {
    const n = heightDiffs.length;
    const measuredEnd = startHeight + heightDiffs.reduce((a, b) => a + b, 0);
    const closure = measuredEnd - endHeight;
    const totalDist = distances.reduce((a, b) => a + b, 0);
    
    const heights: number[] = [startHeight];
    for (let i = 0; i < n; i++) {
      const correction = -closure * distances[i] / totalDist;
      const correctedDiff = heightDiffs[i] + correction;
      heights.push(heights[heights.length - 1] + correctedDiff);
    }
    
    return { heights, closure };
  }

  // ==================== 误差计算 ====================
  
  // 平均误差
  meanError(residuals: number[]): number {
    const n = residuals.length;
    return Math.sqrt(residuals.reduce((sum, r) => sum + r * r, 0) / n);
  }
  
  // 中误差
  standardError(residuals: number[], freedomDegree?: number): number {
    const n = residuals.length;
    const df = freedomDegree || n - 1;
    return Math.sqrt(residuals.reduce((sum, r) => sum + r * r, 0) / df);
  }
  
  // 限差检验
  toleranceCheck(error: number, tolerance: number): boolean {
    return Math.abs(error) <= tolerance;
  }

  // ==================== 辅助功能 ====================
  
  // 3D距离
  distance3D(p1: Point, p2: Point): number {
    const dx = p2.x - p1.x;
    const dy = p2.y - p1.y;
    const dz = (p2.z || 0) - (p1.z || 0);
    return Math.sqrt(dx*dx + dy*dy + dz*dz);
  }
  
  // 斜距转平距
  slopeToHorizontal(slopeDist: number, verticalAngle: number): number {
    return slopeDist * Math.cos(this.degToRad(verticalAngle));
  }
  
  // 三角高程
  trigHeight(distance: number, verticalAngle: number, instrHeight: number, targetHeight: number): number {
    return distance * Math.tan(this.degToRad(verticalAngle)) + instrHeight - targetHeight;
  }
}

export const surveyCalc = new SurveyCalculator();
export default SurveyCalculator;
