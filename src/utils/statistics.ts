/**
 * 统计计算引擎
 */

export interface StatsSummary {
  count: number;
  sum: number;
  mean: number;
  median: number;
  mode: number[];
  min: number;
  max: number;
  range: number;
  variance: number;
  stdDev: number;
  sampleVariance: number;
  sampleStdDev: number;
  skewness: number;
  kurtosis: number;
  q1: number;
  q2: number;
  q3: number;
  iqr: number;
}

export interface RegressionResult {
  slope: number;
  intercept: number;
  r: number;
  rSquared: number;
  equation: string;
  standardError: number;
}

class StatisticsCalculator {
  private dataX: number[] = [];
  private dataY: number[] = [];

  clearData(): void {
    this.dataX = [];
    this.dataY = [];
  }

  setData(data: number[]): void {
    this.dataX = [...data];
    this.dataY = [];
  }

  setDataXY(x: number[], y: number[]): void {
    if (x.length !== y.length) throw new Error('X和Y长度必须相同');
    this.dataX = [...x];
    this.dataY = [...y];
  }

  addData(x: number, y?: number): void {
    this.dataX.push(x);
    if (y !== undefined) this.dataY.push(y);
  }

  getData(): number[] { return [...this.dataX]; }

  // ==================== 基本统计量 ====================
  
  count(): number { return this.dataX.length; }
  
  sum(): number { return this.dataX.reduce((a, b) => a + b, 0); }
  
  sumSquares(): number { return this.dataX.reduce((a, b) => a + b * b, 0); }
  
  mean(): number {
    if (this.dataX.length === 0) return 0;
    return this.sum() / this.dataX.length;
  }
  
  geometricMean(): number {
    if (this.dataX.some(x => x <= 0)) throw new Error('所有数据必须为正');
    const product = this.dataX.reduce((a, b) => a * b, 1);
    return Math.pow(product, 1 / this.dataX.length);
  }
  
  harmonicMean(): number {
    if (this.dataX.some(x => x === 0)) throw new Error('数据不能为零');
    return this.dataX.length / this.dataX.reduce((a, b) => a + 1/b, 0);
  }
  
  median(): number {
    if (this.dataX.length === 0) return 0;
    const sorted = [...this.dataX].sort((a, b) => a - b);
    const mid = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
  }
  
  mode(): number[] {
    if (this.dataX.length === 0) return [];
    const counts = new Map<number, number>();
    let maxCount = 0;
    
    for (const x of this.dataX) {
      const count = (counts.get(x) || 0) + 1;
      counts.set(x, count);
      maxCount = Math.max(maxCount, count);
    }
    
    return Array.from(counts.entries()).filter(([_, c]) => c === maxCount).map(([v]) => v);
  }
  
  min(): number { return Math.min(...this.dataX); }
  max(): number { return Math.max(...this.dataX); }
  range(): number { return this.max() - this.min(); }

  // ==================== 离散度 ====================
  
  variance(): number {
    if (this.dataX.length < 1) return 0;
    const m = this.mean();
    return this.dataX.reduce((sum, x) => sum + (x - m) ** 2, 0) / this.dataX.length;
  }
  
  sampleVariance(): number {
    if (this.dataX.length < 2) return 0;
    const m = this.mean();
    return this.dataX.reduce((sum, x) => sum + (x - m) ** 2, 0) / (this.dataX.length - 1);
  }
  
  stdDev(): number { return Math.sqrt(this.variance()); }
  sampleStdDev(): number { return Math.sqrt(this.sampleVariance()); }
  
  coefficientOfVariation(): number {
    const m = this.mean();
    return m === 0 ? 0 : this.sampleStdDev() / Math.abs(m) * 100;
  }
  
  standardError(): number {
    return this.sampleStdDev() / Math.sqrt(this.dataX.length);
  }
  
  skewness(): number {
    const n = this.dataX.length;
    if (n < 3) return 0;
    const m = this.mean();
    const s = this.sampleStdDev();
    if (s === 0) return 0;
    return (n / ((n-1) * (n-2))) * this.dataX.reduce((sum, x) => sum + ((x-m)/s)**3, 0);
  }
  
  kurtosis(): number {
    const n = this.dataX.length;
    if (n < 4) return 0;
    const m = this.mean();
    const s = this.sampleStdDev();
    if (s === 0) return 0;
    const m4 = this.dataX.reduce((sum, x) => sum + (x-m)**4, 0) / n;
    return m4 / (s ** 4) - 3;
  }

  // ==================== 分位数 ====================
  
  percentile(p: number): number {
    if (this.dataX.length === 0) return 0;
    if (p < 0 || p > 100) throw new Error('百分位数必须在0-100之间');
    
    const sorted = [...this.dataX].sort((a, b) => a - b);
    const k = (sorted.length - 1) * p / 100;
    const f = Math.floor(k);
    const c = Math.ceil(k);
    
    if (f === c) return sorted[f];
    return sorted[f] * (c - k) + sorted[c] * (k - f);
  }
  
  quartiles(): { q1: number; q2: number; q3: number } {
    return {
      q1: this.percentile(25),
      q2: this.percentile(50),
      q3: this.percentile(75)
    };
  }
  
  iqr(): number {
    const { q1, q3 } = this.quartiles();
    return q3 - q1;
  }

  // ==================== 回归分析 ====================
  
  linearRegression(): RegressionResult {
    const n = this.dataX.length;
    if (n !== this.dataY.length || n < 2) throw new Error('需要至少2对数据');
    
    const sumX = this.dataX.reduce((a, b) => a + b, 0);
    const sumY = this.dataY.reduce((a, b) => a + b, 0);
    const sumXY = this.dataX.reduce((sum, x, i) => sum + x * this.dataY[i], 0);
    const sumX2 = this.dataX.reduce((sum, x) => sum + x * x, 0);
    const sumY2 = this.dataY.reduce((sum, y) => sum + y * y, 0);
    
    const meanX = sumX / n;
    const meanY = sumY / n;
    
    const denom = n * sumX2 - sumX * sumX;
    if (Math.abs(denom) < 1e-15) throw new Error('无法计算回归系数');
    
    const slope = (n * sumXY - sumX * sumY) / denom;
    const intercept = meanY - slope * meanX;
    
    const denomR = Math.sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
    const r = Math.abs(denomR) > 1e-15 ? (n * sumXY - sumX * sumY) / denomR : 0;
    
    // 预测值和残差
    const predicted = this.dataX.map(x => intercept + slope * x);
    const sse = this.dataY.reduce((sum, y, i) => sum + (y - predicted[i]) ** 2, 0);
    const se = n > 2 ? Math.sqrt(sse / (n - 2)) : 0;
    
    return {
      slope,
      intercept,
      r,
      rSquared: r * r,
      equation: `y = ${intercept.toFixed(6)} + ${slope.toFixed(6)}x`,
      standardError: se
    };
  }
  
  exponentialRegression(): RegressionResult {
    if (this.dataY.some(y => y <= 0)) throw new Error('Y值必须为正');
    
    const logY = this.dataY.map(y => Math.log(y));
    const origY = this.dataY;
    this.dataY = logY;
    
    const result = this.linearRegression();
    this.dataY = origY;
    
    const a = Math.exp(result.intercept);
    const b = result.slope;
    
    return {
      slope: b,
      intercept: a,
      r: result.r,
      rSquared: result.rSquared,
      equation: `y = ${a.toFixed(6)} × e^(${b.toFixed(6)}x)`,
      standardError: result.standardError
    };
  }
  
  logarithmicRegression(): RegressionResult {
    if (this.dataX.some(x => x <= 0)) throw new Error('X值必须为正');
    
    const logX = this.dataX.map(x => Math.log(x));
    const origX = this.dataX;
    this.dataX = logX;
    
    const result = this.linearRegression();
    this.dataX = origX;
    
    return {
      slope: result.slope,
      intercept: result.intercept,
      r: result.r,
      rSquared: result.rSquared,
      equation: `y = ${result.intercept.toFixed(6)} + ${result.slope.toFixed(6)} × ln(x)`,
      standardError: result.standardError
    };
  }
  
  powerRegression(): RegressionResult {
    if (this.dataX.some(x => x <= 0) || this.dataY.some(y => y <= 0)) 
      throw new Error('X和Y值必须为正');
    
    const logX = this.dataX.map(x => Math.log(x));
    const logY = this.dataY.map(y => Math.log(y));
    const origX = this.dataX;
    const origY = this.dataY;
    this.dataX = logX;
    this.dataY = logY;
    
    const result = this.linearRegression();
    this.dataX = origX;
    this.dataY = origY;
    
    const a = Math.exp(result.intercept);
    const b = result.slope;
    
    return {
      slope: b,
      intercept: a,
      r: result.r,
      rSquared: result.rSquared,
      equation: `y = ${a.toFixed(6)} × x^${b.toFixed(6)}`,
      standardError: result.standardError
    };
  }
  
  correlation(): number {
    return this.linearRegression().r;
  }

  // ==================== 概率分布 ====================
  
  normalPdf(x: number, mean: number = 0, std: number = 1): number {
    return Math.exp(-0.5 * ((x - mean) / std) ** 2) / (std * Math.sqrt(2 * Math.PI));
  }
  
  normalCdf(x: number, mean: number = 0, std: number = 1): number {
    const z = (x - mean) / std;
    // 误差函数近似
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;
    
    const sign = z >= 0 ? 1 : -1;
    const t = 1.0 / (1.0 + p * Math.abs(z));
    const y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.exp(-z*z/2);
    
    return 0.5 * (1 + sign * y);
  }
  
  binomialPmf(k: number, n: number, p: number): number {
    if (k < 0 || k > n) return 0;
    const comb = this.combination(n, k);
    return comb * Math.pow(p, k) * Math.pow(1 - p, n - k);
  }
  
  private combination(n: number, k: number): number {
    if (k > n) return 0;
    if (k === 0 || k === n) return 1;
    let result = 1;
    for (let i = 0; i < k; i++) {
      result *= (n - i) / (i + 1);
    }
    return Math.round(result);
  }
  
  poissonPmf(k: number, lambda: number): number {
    if (k < 0) return 0;
    return Math.pow(lambda, k) * Math.exp(-lambda) / this.factorial(k);
  }
  
  private factorial(n: number): number {
    if (n <= 1) return 1;
    let result = 1;
    for (let i = 2; i <= n; i++) result *= i;
    return result;
  }

  // ==================== 统计摘要 ====================
  
  summary(): StatsSummary {
    const { q1, q2, q3 } = this.quartiles();
    return {
      count: this.count(),
      sum: this.sum(),
      mean: this.mean(),
      median: this.median(),
      mode: this.mode(),
      min: this.min(),
      max: this.max(),
      range: this.range(),
      variance: this.variance(),
      stdDev: this.stdDev(),
      sampleVariance: this.sampleVariance(),
      sampleStdDev: this.sampleStdDev(),
      skewness: this.skewness(),
      kurtosis: this.kurtosis(),
      q1, q2, q3,
      iqr: this.iqr()
    };
  }
}

export const statsCalc = new StatisticsCalculator();
export default StatisticsCalculator;
