/**
 * 科学计算引擎
 * 提供完整的科学计算功能
 */

// 角度单位类型
export type AngleUnit = 'DEG' | 'RAD' | 'GRAD';

// 计算结果类型
export interface CalcResult {
  value: number | string | number[];
  display: string;
  error?: string;
}

// 历史记录类型
export interface HistoryItem {
  id: string;
  expression: string;
  result: string;
  timestamp: number;
  type: 'basic' | 'scientific' | 'survey' | 'statistics' | 'matrix';
}

// 变量存储
export interface Variable {
  name: string;
  value: number;
  description?: string;
}

class ScientificCalculator {
  private angleUnit: AngleUnit = 'DEG';
  private precision: number = 10;
  private memory: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // M1-M10
  private variables: Map<string, number> = new Map();
  private ans: number = 0;

  // 数学常数
  static readonly CONSTANTS = {
    PI: Math.PI,
    E: Math.E,
    PHI: (1 + Math.sqrt(5)) / 2, // 黄金比例
    C: 299792458, // 光速 m/s
    G: 9.80665, // 重力加速度
    H: 6.62607015e-34, // 普朗克常数
    NA: 6.02214076e23, // 阿伏伽德罗常数
    R: 8.31446261815324, // 气体常数
    KB: 1.380649e-23, // 玻尔兹曼常数
  };

  setAngleUnit(unit: AngleUnit): void {
    this.angleUnit = unit;
  }

  getAngleUnit(): AngleUnit {
    return this.angleUnit;
  }

  setPrecision(precision: number): void {
    this.precision = Math.min(Math.max(precision, 0), 15);
  }

  // 角度转换
  toRadians(angle: number): number {
    switch (this.angleUnit) {
      case 'DEG': return angle * Math.PI / 180;
      case 'GRAD': return angle * Math.PI / 200;
      default: return angle;
    }
  }

  fromRadians(angle: number): number {
    switch (this.angleUnit) {
      case 'DEG': return angle * 180 / Math.PI;
      case 'GRAD': return angle * 200 / Math.PI;
      default: return angle;
    }
  }

  // 格式化结果
  formatResult(value: number): string {
    if (!isFinite(value)) {
      return isNaN(value) ? 'Error' : (value > 0 ? '∞' : '-∞');
    }
    if (Math.abs(value) < 1e-10) return '0';
    if (Math.abs(value) >= 1e10 || (Math.abs(value) < 1e-6 && value !== 0)) {
      return value.toExponential(this.precision);
    }
    return parseFloat(value.toPrecision(this.precision)).toString();
  }

  // ==================== 基本运算 ====================
  
  add(a: number, b: number): number { return a + b; }
  subtract(a: number, b: number): number { return a - b; }
  multiply(a: number, b: number): number { return a * b; }
  divide(a: number, b: number): number {
    if (b === 0) throw new Error('除数不能为零');
    return a / b;
  }
  modulo(a: number, b: number): number { return a % b; }
  power(base: number, exp: number): number { return Math.pow(base, exp); }
  sqrt(x: number): number {
    if (x < 0) throw new Error('负数不能开平方根');
    return Math.sqrt(x);
  }
  cbrt(x: number): number { return Math.cbrt(x); }
  nthRoot(x: number, n: number): number {
    if (n === 0) throw new Error('n不能为零');
    if (x < 0 && n % 2 === 0) throw new Error('负数不能开偶数次方根');
    return x >= 0 ? Math.pow(x, 1/n) : -Math.pow(-x, 1/n);
  }
  reciprocal(x: number): number {
    if (x === 0) throw new Error('零没有倒数');
    return 1 / x;
  }
  abs(x: number): number { return Math.abs(x); }
  negate(x: number): number { return -x; }
  square(x: number): number { return x * x; }
  cube(x: number): number { return x * x * x; }

  // ==================== 三角函数 ====================
  
  sin(x: number): number { return Math.sin(this.toRadians(x)); }
  cos(x: number): number { return Math.cos(this.toRadians(x)); }
  tan(x: number): number {
    const rad = this.toRadians(x);
    if (Math.abs(Math.cos(rad)) < 1e-15) throw new Error('正切无定义');
    return Math.tan(rad);
  }
  asin(x: number): number {
    if (Math.abs(x) > 1) throw new Error('反正弦参数必须在[-1,1]');
    return this.fromRadians(Math.asin(x));
  }
  acos(x: number): number {
    if (Math.abs(x) > 1) throw new Error('反余弦参数必须在[-1,1]');
    return this.fromRadians(Math.acos(x));
  }
  atan(x: number): number { return this.fromRadians(Math.atan(x)); }
  atan2(y: number, x: number): number { return this.fromRadians(Math.atan2(y, x)); }
  
  sec(x: number): number { return 1 / this.cos(x); }
  csc(x: number): number { return 1 / this.sin(x); }
  cot(x: number): number { return 1 / this.tan(x); }
  
  // 反三角函数
  asec(x: number): number {
    if (Math.abs(x) < 1) throw new Error('参数绝对值必须≥1');
    return this.acos(1/x);
  }
  acsc(x: number): number {
    if (Math.abs(x) < 1) throw new Error('参数绝对值必须≥1');
    return this.asin(1/x);
  }
  acot(x: number): number { return this.atan(1/x); }

  // ==================== 双曲函数 ====================
  
  sinh(x: number): number { return Math.sinh(x); }
  cosh(x: number): number { return Math.cosh(x); }
  tanh(x: number): number { return Math.tanh(x); }
  asinh(x: number): number { return Math.asinh(x); }
  acosh(x: number): number {
    if (x < 1) throw new Error('反双曲余弦参数必须≥1');
    return Math.acosh(x);
  }
  atanh(x: number): number {
    if (Math.abs(x) >= 1) throw new Error('反双曲正切参数必须在(-1,1)');
    return Math.atanh(x);
  }
  sech(x: number): number { return 1 / Math.cosh(x); }
  csch(x: number): number { return 1 / Math.sinh(x); }
  coth(x: number): number { return 1 / Math.tanh(x); }

  // ==================== 对数和指数 ====================
  
  ln(x: number): number {
    if (x <= 0) throw new Error('对数参数必须>0');
    return Math.log(x);
  }
  log10(x: number): number {
    if (x <= 0) throw new Error('对数参数必须>0');
    return Math.log10(x);
  }
  log2(x: number): number {
    if (x <= 0) throw new Error('对数参数必须>0');
    return Math.log2(x);
  }
  log(x: number, base: number): number {
    if (x <= 0 || base <= 0 || base === 1) throw new Error('无效的对数参数');
    return Math.log(x) / Math.log(base);
  }
  exp(x: number): number { return Math.exp(x); }
  exp10(x: number): number { return Math.pow(10, x); }
  exp2(x: number): number { return Math.pow(2, x); }

  // ==================== 阶乘和组合 ====================
  
  factorial(n: number): number {
    if (n < 0 || !Number.isInteger(n)) throw new Error('阶乘参数必须是非负整数');
    if (n > 170) throw new Error('数值太大');
    let result = 1;
    for (let i = 2; i <= n; i++) result *= i;
    return result;
  }

  permutation(n: number, r: number): number {
    if (n < 0 || r < 0 || !Number.isInteger(n) || !Number.isInteger(r)) 
      throw new Error('参数必须是非负整数');
    if (r > n) throw new Error('r不能大于n');
    return this.factorial(n) / this.factorial(n - r);
  }

  combination(n: number, r: number): number {
    if (n < 0 || r < 0 || !Number.isInteger(n) || !Number.isInteger(r)) 
      throw new Error('参数必须是非负整数');
    if (r > n) throw new Error('r不能大于n');
    return this.factorial(n) / (this.factorial(r) * this.factorial(n - r));
  }

  // ==================== 特殊函数 ====================
  
  gamma(x: number): number {
    // Lanczos近似
    if (x <= 0 && Number.isInteger(x)) throw new Error('负整数没有伽马函数值');
    
    const g = 7;
    const c = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7
    ];
    
    if (x < 0.5) {
      return Math.PI / (Math.sin(Math.PI * x) * this.gamma(1 - x));
    }
    
    x -= 1;
    let a = c[0];
    const t = x + g + 0.5;
    
    for (let i = 1; i < g + 2; i++) {
      a += c[i] / (x + i);
    }
    
    return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
  }

  erf(x: number): number {
    // 近似计算误差函数
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;

    const sign = x >= 0 ? 1 : -1;
    x = Math.abs(x);
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y;
  }

  erfc(x: number): number { return 1 - this.erf(x); }

  // ==================== 数值积分 ====================
  
  integrate(f: (x: number) => number, a: number, b: number, n: number = 1000): number {
    // 辛普森法则
    if (n % 2 !== 0) n++;
    const h = (b - a) / n;
    let sum = f(a) + f(b);
    
    for (let i = 1; i < n; i++) {
      const x = a + i * h;
      sum += (i % 2 === 0 ? 2 : 4) * f(x);
    }
    
    return sum * h / 3;
  }

  // 数值微分
  derivative(f: (x: number) => number, x: number, h: number = 1e-8): number {
    return (f(x + h) - f(x - h)) / (2 * h);
  }

  secondDerivative(f: (x: number) => number, x: number, h: number = 1e-5): number {
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
  }

  // ==================== 方程求解 ====================
  
  solveQuadratic(a: number, b: number, c: number): { x1: number | { re: number; im: number }; x2: number | { re: number; im: number } } {
    if (a === 0) throw new Error('a不能为零');
    
    const discriminant = b * b - 4 * a * c;
    
    if (discriminant >= 0) {
      const sqrtD = Math.sqrt(discriminant);
      return {
        x1: (-b + sqrtD) / (2 * a),
        x2: (-b - sqrtD) / (2 * a)
      };
    } else {
      const sqrtD = Math.sqrt(-discriminant);
      return {
        x1: { re: -b / (2 * a), im: sqrtD / (2 * a) },
        x2: { re: -b / (2 * a), im: -sqrtD / (2 * a) }
      };
    }
  }

  // 牛顿法求根
  newtonRaphson(f: (x: number) => number, x0: number, tol: number = 1e-10, maxIter: number = 100): number {
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      const fx = f(x);
      if (Math.abs(fx) < tol) return x;
      const dfx = this.derivative(f, x);
      if (Math.abs(dfx) < 1e-15) throw new Error('导数接近零');
      const xNew = x - fx / dfx;
      if (Math.abs(xNew - x) < tol) return xNew;
      x = xNew;
    }
    throw new Error('未收敛');
  }

  // ==================== 进制转换 ====================
  
  toBinary(n: number): string {
    if (!Number.isInteger(n)) throw new Error('必须是整数');
    return (n >= 0 ? '' : '-') + Math.abs(n).toString(2);
  }
  toOctal(n: number): string {
    if (!Number.isInteger(n)) throw new Error('必须是整数');
    return (n >= 0 ? '' : '-') + Math.abs(n).toString(8);
  }
  toHex(n: number): string {
    if (!Number.isInteger(n)) throw new Error('必须是整数');
    return (n >= 0 ? '' : '-') + Math.abs(n).toString(16).toUpperCase();
  }
  fromBinary(s: string): number { return parseInt(s, 2); }
  fromOctal(s: string): number { return parseInt(s, 8); }
  fromHex(s: string): number { return parseInt(s, 16); }

  // ==================== 位运算 ====================
  
  bitAnd(a: number, b: number): number { return a & b; }
  bitOr(a: number, b: number): number { return a | b; }
  bitXor(a: number, b: number): number { return a ^ b; }
  bitNot(a: number): number { return ~a; }
  leftShift(a: number, n: number): number { return a << n; }
  rightShift(a: number, n: number): number { return a >> n; }

  // ==================== 内存操作 ====================
  
  memoryStore(index: number, value: number): void {
    if (index >= 0 && index < 10) this.memory[index] = value;
  }
  memoryRecall(index: number): number {
    return (index >= 0 && index < 10) ? this.memory[index] : 0;
  }
  memoryAdd(index: number, value: number): void {
    if (index >= 0 && index < 10) this.memory[index] += value;
  }
  memorySubtract(index: number, value: number): void {
    if (index >= 0 && index < 10) this.memory[index] -= value;
  }
  memoryClear(index: number): void {
    if (index >= 0 && index < 10) this.memory[index] = 0;
  }
  memoryClearAll(): void {
    this.memory.fill(0);
  }

  // 变量操作
  setVariable(name: string, value: number): void {
    this.variables.set(name.toUpperCase(), value);
  }
  getVariable(name: string): number {
    return this.variables.get(name.toUpperCase()) || 0;
  }
  getAllVariables(): Map<string, number> {
    return new Map(this.variables);
  }
  clearVariables(): void {
    this.variables.clear();
  }

  // ANS操作
  setAns(value: number): void { this.ans = value; }
  getAns(): number { return this.ans; }

  // ==================== 单位换算 ====================
  
  convertLength(value: number, from: string, to: string): number {
    const toMeter: { [key: string]: number } = {
      'm': 1, 'km': 1000, 'cm': 0.01, 'mm': 0.001,
      'in': 0.0254, 'ft': 0.3048, 'yd': 0.9144, 'mi': 1609.344,
      'nm': 1852, 'um': 1e-6
    };
    if (!(from in toMeter) || !(to in toMeter)) throw new Error('不支持的单位');
    return value * toMeter[from] / toMeter[to];
  }

  convertAngle(value: number, from: string, to: string): number {
    const toRad: { [key: string]: number } = {
      'rad': 1, 'deg': Math.PI / 180, 'grad': Math.PI / 200,
      'turn': 2 * Math.PI, 'arcmin': Math.PI / 10800, 'arcsec': Math.PI / 648000
    };
    if (!(from in toRad) || !(to in toRad)) throw new Error('不支持的单位');
    return value * toRad[from] / toRad[to];
  }

  convertArea(value: number, from: string, to: string): number {
    const toSqm: { [key: string]: number } = {
      'm2': 1, 'km2': 1e6, 'cm2': 1e-4, 'mm2': 1e-6,
      'ha': 1e4, 'mu': 666.667, 'acre': 4046.86, 'sqft': 0.092903
    };
    if (!(from in toSqm) || !(to in toSqm)) throw new Error('不支持的单位');
    return value * toSqm[from] / toSqm[to];
  }

  convertTemperature(value: number, from: string, to: string): number {
    let kelvin: number;
    switch (from) {
      case 'K': kelvin = value; break;
      case 'C': kelvin = value + 273.15; break;
      case 'F': kelvin = (value - 32) * 5 / 9 + 273.15; break;
      default: throw new Error('不支持的单位');
    }
    switch (to) {
      case 'K': return kelvin;
      case 'C': return kelvin - 273.15;
      case 'F': return (kelvin - 273.15) * 9 / 5 + 32;
      default: throw new Error('不支持的单位');
    }
  }
}

export const calculator = new ScientificCalculator();
export default ScientificCalculator;
