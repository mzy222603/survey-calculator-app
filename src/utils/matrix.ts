/**
 * 矩阵计算引擎
 */

export type Matrix = number[][];

class MatrixCalculator {
  private matrices: Map<string, Matrix> = new Map();

  // 存储矩阵
  store(name: string, matrix: Matrix): void {
    this.matrices.set(name.toUpperCase(), matrix.map(row => [...row]));
  }

  get(name: string): Matrix | undefined {
    const m = this.matrices.get(name.toUpperCase());
    return m ? m.map(row => [...row]) : undefined;
  }

  clear(): void {
    this.matrices.clear();
  }

  // 创建矩阵
  create(rows: number, cols: number, values?: number[][]): Matrix {
    if (values) {
      if (values.length !== rows || values.some(r => r.length !== cols)) {
        throw new Error('维度不匹配');
      }
      return values.map(row => [...row]);
    }
    return Array(rows).fill(0).map(() => Array(cols).fill(0));
  }

  identity(n: number): Matrix {
    return Array(n).fill(0).map((_, i) => Array(n).fill(0).map((_, j) => i === j ? 1 : 0));
  }

  zeros(rows: number, cols: number): Matrix {
    return this.create(rows, cols);
  }

  ones(rows: number, cols: number): Matrix {
    return Array(rows).fill(0).map(() => Array(cols).fill(1));
  }

  diagonal(values: number[]): Matrix {
    const n = values.length;
    return Array(n).fill(0).map((_, i) => Array(n).fill(0).map((_, j) => i === j ? values[i] : 0));
  }

  // 获取维度
  dimensions(m: Matrix): { rows: number; cols: number } {
    return { rows: m.length, cols: m[0]?.length || 0 };
  }

  // ==================== 基本运算 ====================

  add(a: Matrix, b: Matrix): Matrix {
    const { rows, cols } = this.dimensions(a);
    const dimB = this.dimensions(b);
    if (rows !== dimB.rows || cols !== dimB.cols) throw new Error('维度不匹配');
    return a.map((row, i) => row.map((val, j) => val + b[i][j]));
  }

  subtract(a: Matrix, b: Matrix): Matrix {
    const { rows, cols } = this.dimensions(a);
    const dimB = this.dimensions(b);
    if (rows !== dimB.rows || cols !== dimB.cols) throw new Error('维度不匹配');
    return a.map((row, i) => row.map((val, j) => val - b[i][j]));
  }

  multiply(a: Matrix, b: Matrix): Matrix {
    const { rows: rowsA, cols: colsA } = this.dimensions(a);
    const { rows: rowsB, cols: colsB } = this.dimensions(b);
    if (colsA !== rowsB) throw new Error(`维度不匹配: ${rowsA}×${colsA} × ${rowsB}×${colsB}`);
    
    const result = this.zeros(rowsA, colsB);
    for (let i = 0; i < rowsA; i++) {
      for (let j = 0; j < colsB; j++) {
        for (let k = 0; k < colsA; k++) {
          result[i][j] += a[i][k] * b[k][j];
        }
      }
    }
    return result;
  }

  scalarMultiply(m: Matrix, scalar: number): Matrix {
    return m.map(row => row.map(val => val * scalar));
  }

  transpose(m: Matrix): Matrix {
    const { rows, cols } = this.dimensions(m);
    return Array(cols).fill(0).map((_, i) => Array(rows).fill(0).map((_, j) => m[j][i]));
  }

  // ==================== 行列式和逆矩阵 ====================

  determinant(m: Matrix): number {
    const { rows, cols } = this.dimensions(m);
    if (rows !== cols) throw new Error('只有方阵有行列式');
    
    const n = rows;
    if (n === 1) return m[0][0];
    if (n === 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    
    // LU分解法
    const mat = m.map(row => [...row]);
    let det = 1;
    
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(mat[k][i]) > Math.abs(mat[maxRow][i])) maxRow = k;
      }
      
      if (maxRow !== i) {
        [mat[i], mat[maxRow]] = [mat[maxRow], mat[i]];
        det *= -1;
      }
      
      if (Math.abs(mat[i][i]) < 1e-15) return 0;
      det *= mat[i][i];
      
      for (let k = i + 1; k < n; k++) {
        const factor = mat[k][i] / mat[i][i];
        for (let j = i; j < n; j++) {
          mat[k][j] -= factor * mat[i][j];
        }
      }
    }
    
    return det;
  }

  inverse(m: Matrix): Matrix {
    const { rows, cols } = this.dimensions(m);
    if (rows !== cols) throw new Error('只有方阵可逆');
    
    const n = rows;
    const augmented = m.map((row, i) => [...row, ...this.identity(n)[i]]);
    
    // 高斯-约旦消元
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) maxRow = k;
      }
      
      [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];
      
      if (Math.abs(augmented[i][i]) < 1e-15) throw new Error('矩阵不可逆');
      
      const pivot = augmented[i][i];
      for (let j = 0; j < 2 * n; j++) augmented[i][j] /= pivot;
      
      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const factor = augmented[k][i];
          for (let j = 0; j < 2 * n; j++) {
            augmented[k][j] -= factor * augmented[i][j];
          }
        }
      }
    }
    
    return augmented.map(row => row.slice(n));
  }

  rank(m: Matrix): number {
    const mat = m.map(row => [...row]);
    const { rows, cols } = this.dimensions(mat);
    let r = 0;
    
    for (let col = 0; col < Math.min(rows, cols); col++) {
      let pivotRow = null;
      for (let row = r; row < rows; row++) {
        if (Math.abs(mat[row][col]) > 1e-15) {
          pivotRow = row;
          break;
        }
      }
      
      if (pivotRow === null) continue;
      
      [mat[r], mat[pivotRow]] = [mat[pivotRow], mat[r]];
      
      for (let row = r + 1; row < rows; row++) {
        if (Math.abs(mat[row][col]) > 1e-15) {
          const factor = mat[row][col] / mat[r][col];
          for (let j = col; j < cols; j++) {
            mat[row][j] -= factor * mat[r][j];
          }
        }
      }
      
      r++;
    }
    
    return r;
  }

  trace(m: Matrix): number {
    const { rows, cols } = this.dimensions(m);
    if (rows !== cols) throw new Error('只有方阵有迹');
    return m.reduce((sum, row, i) => sum + row[i], 0);
  }

  // ==================== 特征值 ====================

  eigenvalues2x2(m: Matrix): { lambda1: number | { re: number; im: number }; lambda2: number | { re: number; im: number } } {
    const { rows, cols } = this.dimensions(m);
    if (rows !== 2 || cols !== 2) throw new Error('只支持2×2矩阵');
    
    const tr = m[0][0] + m[1][1];
    const det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    const disc = tr * tr - 4 * det;
    
    if (disc >= 0) {
      const sqrtDisc = Math.sqrt(disc);
      return {
        lambda1: (tr + sqrtDisc) / 2,
        lambda2: (tr - sqrtDisc) / 2
      };
    } else {
      const sqrtDisc = Math.sqrt(-disc);
      return {
        lambda1: { re: tr / 2, im: sqrtDisc / 2 },
        lambda2: { re: tr / 2, im: -sqrtDisc / 2 }
      };
    }
  }

  // 幂迭代法求最大特征值
  powerIteration(m: Matrix, maxIter: number = 1000, tol: number = 1e-10): { eigenvalue: number; eigenvector: number[] } {
    const n = m.length;
    let v = Array(n).fill(1 / Math.sqrt(n));
    let eigenvalue = 0;
    
    for (let iter = 0; iter < maxIter; iter++) {
      // Av
      const newV = Array(n).fill(0);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          newV[i] += m[i][j] * v[j];
        }
      }
      
      // 归一化
      const norm = Math.sqrt(newV.reduce((s, x) => s + x * x, 0));
      if (norm < 1e-15) throw new Error('迭代失败');
      
      for (let i = 0; i < n; i++) newV[i] /= norm;
      
      // 计算特征值
      const newEigenvalue = newV.reduce((s, _, i) => {
        let sum = 0;
        for (let j = 0; j < n; j++) sum += m[i][j] * v[j];
        return s + newV[i] * sum;
      }, 0);
      
      if (Math.abs(newEigenvalue - eigenvalue) < tol) {
        return { eigenvalue: newEigenvalue, eigenvector: newV };
      }
      
      eigenvalue = newEigenvalue;
      v = newV;
    }
    
    return { eigenvalue, eigenvector: v };
  }

  // ==================== 矩阵分解 ====================

  luDecomposition(m: Matrix): { L: Matrix; U: Matrix } {
    const n = m.length;
    const L = this.identity(n);
    const U = this.zeros(n, n);
    
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let sum = 0;
        for (let k = 0; k < i; k++) sum += L[i][k] * U[k][j];
        U[i][j] = m[i][j] - sum;
      }
      
      for (let j = i + 1; j < n; j++) {
        let sum = 0;
        for (let k = 0; k < i; k++) sum += L[j][k] * U[k][i];
        if (Math.abs(U[i][i]) < 1e-15) throw new Error('LU分解失败');
        L[j][i] = (m[j][i] - sum) / U[i][i];
      }
    }
    
    return { L, U };
  }

  cholesky(m: Matrix): Matrix {
    const n = m.length;
    const L = this.zeros(n, n);
    
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let k = 0; k < j; k++) sum += L[i][k] * L[j][k];
        
        if (i === j) {
          const val = m[i][i] - sum;
          if (val < 0) throw new Error('矩阵不是正定的');
          L[i][j] = Math.sqrt(val);
        } else {
          L[i][j] = (m[i][j] - sum) / L[j][j];
        }
      }
    }
    
    return L;
  }

  // ==================== 线性方程组 ====================

  solveLinearSystem(a: Matrix, b: number[]): number[] {
    const n = a.length;
    const augmented = a.map((row, i) => [...row, b[i]]);
    
    // 前向消元
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) maxRow = k;
      }
      
      [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];
      
      if (Math.abs(augmented[i][i]) < 1e-15) throw new Error('方程组无唯一解');
      
      for (let k = i + 1; k < n; k++) {
        const factor = augmented[k][i] / augmented[i][i];
        for (let j = i; j <= n; j++) {
          augmented[k][j] -= factor * augmented[i][j];
        }
      }
    }
    
    // 回代
    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = augmented[i][n];
      for (let j = i + 1; j < n; j++) {
        x[i] -= augmented[i][j] * x[j];
      }
      x[i] /= augmented[i][i];
    }
    
    return x;
  }

  leastSquares(a: Matrix, b: number[]): number[] {
    const at = this.transpose(a);
    const ata = this.multiply(at, a);
    const atb = at.map(row => row.reduce((s, v, i) => s + v * b[i], 0));
    return this.solveLinearSystem(ata, atb);
  }

  // ==================== 范数 ====================

  normFrobenius(m: Matrix): number {
    return Math.sqrt(m.reduce((s, row) => s + row.reduce((rs, v) => rs + v * v, 0), 0));
  }

  norm1(m: Matrix): number {
    const { cols } = this.dimensions(m);
    let max = 0;
    for (let j = 0; j < cols; j++) {
      const colSum = m.reduce((s, row) => s + Math.abs(row[j]), 0);
      max = Math.max(max, colSum);
    }
    return max;
  }

  normInf(m: Matrix): number {
    return Math.max(...m.map(row => row.reduce((s, v) => s + Math.abs(v), 0)));
  }

  // 格式化输出
  toString(m: Matrix, precision: number = 4): string {
    return m.map(row => '[' + row.map(v => v.toFixed(precision).padStart(precision + 6)).join(' ') + ']').join('\n');
  }
}

export const matrixCalc = new MatrixCalculator();
export default MatrixCalculator;
