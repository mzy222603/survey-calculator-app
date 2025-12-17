import React, { useState, useEffect, useCallback } from 'react';
import './App.css';
import { calculator } from './utils/calculator';
import { surveyCalc, Point } from './utils/survey';
import { statsCalc } from './utils/statistics';
import { storage, Settings, HistoryItem } from './utils/storage';

type TabType = 'home' | 'calc' | 'survey' | 'stats' | 'settings';

function App() {
  const [activeTab, setActiveTab] = useState<TabType>('home');
  const [settings, setSettings] = useState<Settings>(storage.getSettings());
  const [showScientific, setShowScientific] = useState(true);
  
  // 计算器状态
  const [display, setDisplay] = useState('0');
  const [expression, setExpression] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [memory, setMemory] = useState(0);
  const [hasMemory, setHasMemory] = useState(false);
  
  // 测绘状态
  const [surveyType, setSurveyType] = useState('forward');
  const [surveyInputs, setSurveyInputs] = useState<{[key: string]: string}>({});
  const [surveyResult, setSurveyResult] = useState('');
  
  // 统计状态
  const [statsData, setStatsData] = useState('');
  const [statsResult, setStatsResult] = useState('');

  useEffect(() => {
    setHistory(storage.getHistory());
    calculator.setAngleUnit(settings.angleUnit);
  }, [settings.angleUnit]);

  const vibrate = useCallback(() => {
    if (settings.vibration && navigator.vibrate) {
      navigator.vibrate(10);
    }
  }, [settings.vibration]);

  const saveToHistory = (expr: string, result: string, type: HistoryItem['type'] = 'calc') => {
    storage.addHistory({ expression: expr, result, type });
    setHistory(storage.getHistory());
  };

  // ==================== 计算器逻辑 ====================
  
  const clearAll = () => {
    vibrate();
    setDisplay('0');
    setExpression('');
  };

  const appendToDisplay = (val: string) => {
    vibrate();
    if (display === '0' && val !== '.') {
      setDisplay(val);
    } else if (display === 'Error') {
      setDisplay(val);
    } else {
      setDisplay(display + val);
    }
  };

  const backspace = () => {
    vibrate();
    if (display.length > 1 && display !== 'Error') {
      setDisplay(display.slice(0, -1));
    } else {
      setDisplay('0');
    }
  };

  const toggleSign = () => {
    vibrate();
    if (display !== '0' && display !== 'Error') {
      setDisplay(display.startsWith('-') ? display.slice(1) : '-' + display);
    }
  };

  const calculate = () => {
    vibrate();
    try {
      // 替换显示符号为计算符号
      let expr = display
        .replace(/×/g, '*')
        .replace(/÷/g, '/')
        .replace(/π/g, `(${Math.PI})`)
        .replace(/e(?![x])/g, `(${Math.E})`)
        .replace(/\^/g, '**')
        .replace(/mod/g, '%')
        .replace(/√\(/g, 'Math.sqrt(')
        .replace(/∛\(/g, 'Math.cbrt(')
        .replace(/sin\(/g, `Math.sin(${settings.angleUnit === 'DEG' ? 'Math.PI/180*' : settings.angleUnit === 'GRAD' ? 'Math.PI/200*' : ''}`)
        .replace(/cos\(/g, `Math.cos(${settings.angleUnit === 'DEG' ? 'Math.PI/180*' : settings.angleUnit === 'GRAD' ? 'Math.PI/200*' : ''}`)
        .replace(/tan\(/g, `Math.tan(${settings.angleUnit === 'DEG' ? 'Math.PI/180*' : settings.angleUnit === 'GRAD' ? 'Math.PI/200*' : ''}`)
        .replace(/asin\(/g, `(${settings.angleUnit === 'DEG' ? '180/Math.PI*' : settings.angleUnit === 'GRAD' ? '200/Math.PI*' : ''}Math.asin(`)
        .replace(/acos\(/g, `(${settings.angleUnit === 'DEG' ? '180/Math.PI*' : settings.angleUnit === 'GRAD' ? '200/Math.PI*' : ''}Math.acos(`)
        .replace(/atan\(/g, `(${settings.angleUnit === 'DEG' ? '180/Math.PI*' : settings.angleUnit === 'GRAD' ? '200/Math.PI*' : ''}Math.atan(`)
        .replace(/sinh\(/g, 'Math.sinh(')
        .replace(/cosh\(/g, 'Math.cosh(')
        .replace(/tanh\(/g, 'Math.tanh(')
        .replace(/ln\(/g, 'Math.log(')
        .replace(/log\(/g, 'Math.log10(')
        .replace(/log2\(/g, 'Math.log2(')
        .replace(/abs\(/g, 'Math.abs(')
        .replace(/exp\(/g, 'Math.exp(')
        .replace(/floor\(/g, 'Math.floor(')
        .replace(/ceil\(/g, 'Math.ceil(')
        .replace(/round\(/g, 'Math.round(');
      
      // eslint-disable-next-line no-eval
      const result = eval(expr);
      const formatted = calculator.formatResult(result);
      saveToHistory(display, formatted);
      setExpression(display + ' =');
      setDisplay(formatted);
      calculator.setAns(result);
    } catch (e) {
      setDisplay('Error');
    }
  };

  const applyFunction = (func: string) => {
    vibrate();
    try {
      const value = parseFloat(display);
      if (isNaN(value) && !['π', 'e', 'rand'].includes(func)) {
        setDisplay('Error');
        return;
      }
      
      let result: number;
      switch (func) {
        case 'sin': result = calculator.sin(value); break;
        case 'cos': result = calculator.cos(value); break;
        case 'tan': result = calculator.tan(value); break;
        case 'asin': result = calculator.asin(value); break;
        case 'acos': result = calculator.acos(value); break;
        case 'atan': result = calculator.atan(value); break;
        case 'sinh': result = Math.sinh(value); break;
        case 'cosh': result = Math.cosh(value); break;
        case 'tanh': result = Math.tanh(value); break;
        case 'ln': result = Math.log(value); break;
        case 'log': result = Math.log10(value); break;
        case 'log2': result = Math.log2(value); break;
        case '√': result = Math.sqrt(value); break;
        case '∛': result = Math.cbrt(value); break;
        case 'x²': result = value * value; break;
        case 'x³': result = value * value * value; break;
        case '1/x': result = 1 / value; break;
        case 'n!': result = calculator.factorial(Math.round(value)); break;
        case 'abs': result = Math.abs(value); break;
        case 'exp': result = Math.exp(value); break;
        case '10ˣ': result = Math.pow(10, value); break;
        case '2ˣ': result = Math.pow(2, value); break;
        case 'eˣ': result = Math.exp(value); break;
        case 'π': result = Math.PI; break;
        case 'e': result = Math.E; break;
        case 'rand': result = Math.random(); break;
        case 'floor': result = Math.floor(value); break;
        case 'ceil': result = Math.ceil(value); break;
        case 'round': result = Math.round(value); break;
        case '%': result = value / 100; break;
        case 'ANS': result = calculator.getAns(); break;
        default: return;
      }
      
      saveToHistory(`${func}(${display})`, calculator.formatResult(result));
      setDisplay(calculator.formatResult(result));
    } catch (e) {
      setDisplay('Error');
    }
  };

  const insertFunction = (func: string) => {
    vibrate();
    if (display === '0' || display === 'Error') {
      setDisplay(func + '(');
    } else {
      setDisplay(display + func + '(');
    }
  };

  // 内存功能
  const memClear = () => { vibrate(); setMemory(0); setHasMemory(false); };
  const memRecall = () => { vibrate(); setDisplay(String(memory)); };
  const memAdd = () => { vibrate(); setMemory(memory + parseFloat(display) || 0); setHasMemory(true); };
  const memSub = () => { vibrate(); setMemory(memory - parseFloat(display) || 0); setHasMemory(true); };

  // 度分秒转换功能
  // 十进制度 → 度分秒
  const deg2dms = () => {
    vibrate();
    try {
      const deg = parseFloat(display);
      if (isNaN(deg)) { setDisplay('Error'); return; }
      const sign = deg < 0 ? -1 : 1;
      const absDeg = Math.abs(deg);
      const d = Math.floor(absDeg);
      const mFloat = (absDeg - d) * 60;
      const m = Math.floor(mFloat);
      const s = (mFloat - m) * 60;
      // 显示格式: 30°15'20.1234"
      const result = `${sign < 0 ? '-' : ''}${d}°${m}'${s.toFixed(4)}"`;
      saveToHistory(`D→DMS(${display})`, result);
      setDisplay(result);
    } catch { setDisplay('Error'); }
  };

  // 度分秒 → 十进制度
  const dms2deg = () => {
    vibrate();
    try {
      const input = display.trim();
      let deg = 0;
      
      // 格式1: 30°15'20.5" 或 30°15'20.5 或 30°15' 或 30°
      if (input.includes('°')) {
        const parts = input.replace(/['"′″]/g, ' ').replace('°', ' ').trim().split(/\s+/);
        const sign = input.startsWith('-') ? -1 : 1;
        const d = Math.abs(parseFloat(parts[0])) || 0;
        const m = parseFloat(parts[1]) || 0;
        const s = parseFloat(parts[2]) || 0;
        deg = sign * (d + m / 60 + s / 3600);
      }
      // 格式2: 30.1520 表示 30°15'20" (DD.MMSS格式)
      else if (/^-?\d+\.\d{4,}$/.test(input)) {
        const val = parseFloat(input);
        const sign = val < 0 ? -1 : 1;
        const absVal = Math.abs(val);
        const d = Math.floor(absVal);
        const decimal = absVal - d;
        const mm = Math.floor(decimal * 100);
        const ss = (decimal * 100 - mm) * 100;
        deg = sign * (d + mm / 60 + ss / 3600);
      }
      // 格式3: 纯数字，直接当作度
      else {
        deg = parseFloat(input);
      }
      
      if (isNaN(deg)) { setDisplay('Error'); return; }
      const result = deg.toFixed(8);
      saveToHistory(`DMS→D(${display})`, result);
      setDisplay(result);
    } catch { setDisplay('Error'); }
  };

  // 弧度 → 度
  const rad2deg = () => {
    vibrate();
    try {
      const rad = parseFloat(display);
      if (isNaN(rad)) { setDisplay('Error'); return; }
      const deg = rad * 180 / Math.PI;
      const result = deg.toFixed(8);
      saveToHistory(`R→D(${display})`, result);
      setDisplay(result);
    } catch { setDisplay('Error'); }
  };

  // 度 → 弧度
  const deg2rad = () => {
    vibrate();
    try {
      const deg = parseFloat(display);
      if (isNaN(deg)) { setDisplay('Error'); return; }
      const rad = deg * Math.PI / 180;
      const result = rad.toFixed(8);
      saveToHistory(`D→R(${display})`, result);
      setDisplay(result);
    } catch { setDisplay('Error'); }
  };

  // 插入度分秒符号
  const insertDMS = (symbol: string) => {
    vibrate();
    if (display === '0' || display === 'Error') {
      setDisplay(symbol);
    } else {
      setDisplay(display + symbol);
    }
  };

  // ==================== 测绘计算 ====================
  
  const handleSurveyInput = (key: string, value: string) => {
    setSurveyInputs(prev => ({ ...prev, [key]: value }));
  };

  const calculateSurvey = () => {
    vibrate();
    try {
      let result = '';
      
      switch (surveyType) {
        case 'forward': {
          const x0 = parseFloat(surveyInputs.x0 || '0');
          const y0 = parseFloat(surveyInputs.y0 || '0');
          const azimuth = parseFloat(surveyInputs.azimuth || '0');
          const distance = parseFloat(surveyInputs.distance || '0');
          const p = surveyCalc.forwardCalc({ x: x0, y: y0 }, azimuth, distance);
          result = `【坐标正算结果】

起点坐标:
  X = ${x0.toFixed(4)}
  Y = ${y0.toFixed(4)}

方位角: ${azimuth.toFixed(6)}°
距离: ${distance.toFixed(4)} m

━━━━━━━━━━━━━━
计算点坐标:
  X = ${p.x.toFixed(4)} m
  Y = ${p.y.toFixed(4)} m`;
          break;
        }
        case 'inverse': {
          const x1 = parseFloat(surveyInputs.x1 || '0');
          const y1 = parseFloat(surveyInputs.y1 || '0');
          const x2 = parseFloat(surveyInputs.x2 || '100');
          const y2 = parseFloat(surveyInputs.y2 || '100');
          const inv = surveyCalc.inverseCalc({ x: x1, y: y1 }, { x: x2, y: y2 });
          result = `【坐标反算结果】

起点: (${x1.toFixed(4)}, ${y1.toFixed(4)})
终点: (${x2.toFixed(4)}, ${y2.toFixed(4)})

━━━━━━━━━━━━━━
方位角: ${inv.azimuth.toFixed(6)}°
距离: ${inv.distance.toFixed(4)} m`;
          break;
        }
        case 'forward_intersection': {
          const xa = parseFloat(surveyInputs.xa || '0');
          const ya = parseFloat(surveyInputs.ya || '0');
          const xb = parseFloat(surveyInputs.xb || '100');
          const yb = parseFloat(surveyInputs.yb || '0');
          const angleA = parseFloat(surveyInputs.angleA || '45');
          const angleB = parseFloat(surveyInputs.angleB || '45');
          const p = surveyCalc.forwardIntersection({x: xa, y: ya}, {x: xb, y: yb}, angleA, angleB);
          result = `【前方交会结果】

A点: (${xa.toFixed(4)}, ${ya.toFixed(4)})
B点: (${xb.toFixed(4)}, ${yb.toFixed(4)})
∠PAB: ${angleA.toFixed(6)}°
∠PBA: ${angleB.toFixed(6)}°

━━━━━━━━━━━━━━
P点坐标:
  X = ${p.x.toFixed(4)} m
  Y = ${p.y.toFixed(4)} m`;
          break;
        }
        case 'resection': {
          const xa = parseFloat(surveyInputs.xa || '0');
          const ya = parseFloat(surveyInputs.ya || '0');
          const xb = parseFloat(surveyInputs.xb || '100');
          const yb = parseFloat(surveyInputs.yb || '0');
          const xc = parseFloat(surveyInputs.xc || '50');
          const yc = parseFloat(surveyInputs.yc || '100');
          const alpha = parseFloat(surveyInputs.alpha || '60');
          const beta = parseFloat(surveyInputs.beta || '60');
          const p = surveyCalc.resection({x: xa, y: ya}, {x: xb, y: yb}, {x: xc, y: yc}, alpha, beta);
          result = `【后方交会结果】

A点: (${xa.toFixed(4)}, ${ya.toFixed(4)})
B点: (${xb.toFixed(4)}, ${yb.toFixed(4)})
C点: (${xc.toFixed(4)}, ${yc.toFixed(4)})
∠APB: ${alpha.toFixed(6)}°
∠BPC: ${beta.toFixed(6)}°

━━━━━━━━━━━━━━
测站P坐标:
  X = ${p.x.toFixed(4)} m
  Y = ${p.y.toFixed(4)} m`;
          break;
        }
        case 'side_shot': {
          const x0 = parseFloat(surveyInputs.x0 || '0');
          const y0 = parseFloat(surveyInputs.y0 || '0');
          const backAzimuth = parseFloat(surveyInputs.backAzimuth || '0');
          const angle = parseFloat(surveyInputs.angle || '90');
          const distance = parseFloat(surveyInputs.distance || '100');
          const azimuth = (backAzimuth + angle + 180) % 360;
          const p = surveyCalc.forwardCalc({x: x0, y: y0}, azimuth, distance);
          result = `【侧方交会/支距法】

测站: (${x0.toFixed(4)}, ${y0.toFixed(4)})
后视方位角: ${backAzimuth.toFixed(6)}°
水平角: ${angle.toFixed(6)}°
距离: ${distance.toFixed(4)} m

计算方位角: ${azimuth.toFixed(6)}°

━━━━━━━━━━━━━━
目标点坐标:
  X = ${p.x.toFixed(4)} m
  Y = ${p.y.toFixed(4)} m`;
          break;
        }
        case 'area': {
          const pointsStr = surveyInputs.points || '0,0\n100,0\n100,100\n0,100';
          const points: Point[] = pointsStr.split('\n').map(line => {
            const [x, y] = line.split(',').map(Number);
            return { x: x || 0, y: y || 0 };
          });
          const area = surveyCalc.polygonArea(points);
          let pointsList = points.map((p, i) => `  ${i+1}. (${p.x.toFixed(4)}, ${p.y.toFixed(4)})`).join('\n');
          result = `【多边形面积计算】

顶点坐标:
${pointsList}

顶点数: ${points.length}

━━━━━━━━━━━━━━
面积 = ${area.toFixed(4)} m²
     = ${(area/10000).toFixed(6)} 公顷
     = ${(area/666.67).toFixed(4)} 亩`;
          break;
        }
        case 'gauss_forward': {
          const lat = parseFloat(surveyInputs.lat || '30');
          const lon = parseFloat(surveyInputs.lon || '120');
          const cm = surveyInputs.cm ? parseFloat(surveyInputs.cm) : undefined;
          const g = surveyCalc.gaussForward(lat, lon, cm);
          result = `【高斯投影正算】

经度: ${lon.toFixed(8)}°
纬度: ${lat.toFixed(8)}°
中央子午线: ${g.centralMeridian}°

━━━━━━━━━━━━━━
投影坐标:
  X (北向) = ${g.x.toFixed(4)} m
  Y (东向) = ${g.y.toFixed(4)} m

带号: ${g.zone}`;
          break;
        }
        case 'gauss_inverse': {
          const x = parseFloat(surveyInputs.gx || '3000000');
          const y = parseFloat(surveyInputs.gy || '500000');
          const cm = parseFloat(surveyInputs.cm || '120');
          const g = surveyCalc.gaussInverse(x, y, cm);
          result = `【高斯投影反算】

投影坐标:
  X = ${x.toFixed(4)} m
  Y = ${y.toFixed(4)} m
中央子午线: ${cm}°

━━━━━━━━━━━━━━
地理坐标:
  纬度 = ${g.lat.toFixed(8)}°
  经度 = ${g.lon.toFixed(8)}°`;
          break;
        }
        case 'curve': {
          const radius = parseFloat(surveyInputs.radius || '500');
          const deflection = parseFloat(surveyInputs.deflection || '30');
          const c = surveyCalc.circularCurve(radius, deflection);
          result = `【圆曲线要素计算】

半径 R = ${radius.toFixed(4)} m
转向角 α = ${deflection.toFixed(6)}°

━━━━━━━━━━━━━━
曲线要素:
  切线长 T = ${c.tangentLength.toFixed(4)} m
  曲线长 L = ${c.curveLength.toFixed(4)} m
  外矢距 E = ${c.externalDistance.toFixed(4)} m
  弦长 C = ${c.chord.toFixed(4)} m`;
          break;
        }
        case 'traverse': {
          const startX = parseFloat(surveyInputs.startX || '0');
          const startY = parseFloat(surveyInputs.startY || '0');
          const startAz = parseFloat(surveyInputs.startAz || '0');
          const stationsStr = surveyInputs.stations || '90,100\n90,100\n90,100\n90,100';
          const stations = stationsStr.split('\n').map(line => {
            const [angle, dist] = line.split(',').map(Number);
            return { angle: angle || 0, distance: dist || 0 };
          });
          const tr = surveyCalc.closedTraverse({x: startX, y: startY}, startAz, stations);
          let stationResults = tr.points.map((p, i) =>
            `  ${i+1}. X=${p.x.toFixed(4)}, Y=${p.y.toFixed(4)}`
          ).join('\n');
          result = `【闭合导线计算】

起始点: (${startX}, ${startY})
起始方位角: ${startAz}°
测站数: ${stations.length}

━━━━━━━━━━━━━━
角度闭合差: ${(tr.angleClosure*3600).toFixed(1)}"
相对闭合差: 1/${Math.round(1/tr.relativeClosure)}

平差后坐标:
${stationResults}`;
          break;
        }
        case 'transform4': {
          const dx = parseFloat(surveyInputs.dx || '100');
          const dy = parseFloat(surveyInputs.dy || '200');
          const scale = parseFloat(surveyInputs.scale || '1');
          const rotation = parseFloat(surveyInputs.rotation || '0');
          const x = parseFloat(surveyInputs.tx || '1000');
          const y = parseFloat(surveyInputs.ty || '2000');
          const rot = rotation * Math.PI / 180;
          const newX = dx + scale * (x * Math.cos(rot) - y * Math.sin(rot));
          const newY = dy + scale * (x * Math.sin(rot) + y * Math.cos(rot));
          result = `【四参数坐标转换】

转换参数:
  ΔX = ${dx.toFixed(4)} m
  ΔY = ${dy.toFixed(4)} m
  尺度 = ${scale.toFixed(8)}
  旋转角 = ${rotation.toFixed(8)}°

原坐标: (${x}, ${y})

━━━━━━━━━━━━━━
转换后坐标:
  X' = ${newX.toFixed(4)} m
  Y' = ${newY.toFixed(4)} m`;
          break;
        }
        case 'leveling': {
          const heights = surveyInputs.heights || '100.000\n1.234,-2.345\n1.567,-1.890\n1.123,-2.456';
          const lines = heights.split('\n');
          const startH = parseFloat(lines[0]);
          let h = startH;
          let total = 0;
          let observations: string[] = [];
          for (let i = 1; i < lines.length; i++) {
            const [back, fore] = lines[i].split(',').map(Number);
            const diff = back - fore;
            total += diff;
            h += diff;
            observations.push(`  ${i}. 后视=${back.toFixed(3)}, 前视=${Math.abs(fore).toFixed(3)}, 高差=${diff.toFixed(3)}, H=${h.toFixed(3)}`);
          }
          result = `【水准测量计算】

起始高程: ${startH.toFixed(3)} m
测段数: ${lines.length - 1}

观测数据:
${observations.join('\n')}

━━━━━━━━━━━━━━
总高差: ${total.toFixed(3)} m
终点高程: ${h.toFixed(3)} m`;
          break;
        }
        case 'earthwork': {
          const area1 = parseFloat(surveyInputs.area1 || '100');
          const area2 = parseFloat(surveyInputs.area2 || '120');
          const dist = parseFloat(surveyInputs.edist || '50');
          const avg = (area1 + area2) / 2 * dist;
          const pyramid = dist / 3 * (area1 + area2 + Math.sqrt(area1 * area2));
          result = `【土方量计算】

断面1面积: ${area1.toFixed(4)} m²
断面2面积: ${area2.toFixed(4)} m²
断面间距: ${dist.toFixed(4)} m

━━━━━━━━━━━━━━
平均断面法:
  V = ${avg.toFixed(4)} m³

棱台公式法:
  V = ${pyramid.toFixed(4)} m³`;
          break;
        }
        default:
          result = '请选择计算类型';
      }
      
      setSurveyResult(result);
      saveToHistory(`测绘-${surveyType}`, result.split('\n')[0], 'survey');
    } catch (e: any) {
      setSurveyResult(`计算错误: ${e.message}`);
    }
  };

  // ==================== 统计计算 ====================
  
  const calculateStats = () => {
    vibrate();
    try {
      const values = statsData.split(/[,\s\n]+/).map(Number).filter(n => !isNaN(n));
      if (values.length === 0) {
        setStatsResult('请输入有效数据');
        return;
      }
      
      statsCalc.setData(values);
      const q = statsCalc.quartiles();
      const summary = {
        count: statsCalc.count(),
        sum: statsCalc.sum(),
        mean: statsCalc.mean(),
        median: statsCalc.median(),
        mode: statsCalc.mode(),
        min: statsCalc.min(),
        max: statsCalc.max(),
        range: statsCalc.range(),
        variance: statsCalc.variance(),
        stdDev: statsCalc.stdDev(),
        sampleStdDev: statsCalc.sampleStdDev(),
        skewness: statsCalc.skewness(),
        kurtosis: statsCalc.kurtosis(),
        q1: q.q1,
        q2: q.q2,
        q3: q.q3,
        iqr: statsCalc.iqr()
      };
      
      const result = `【统计分析结果】

样本数: ${summary.count}
━━━━━━━━━━━━━━
集中趋势:
  总和: ${summary.sum.toFixed(6)}
  均值: ${summary.mean.toFixed(6)}
  中位数: ${summary.median.toFixed(6)}
  众数: ${summary.mode.join(', ')}

离散程度:
  最小值: ${summary.min.toFixed(6)}
  最大值: ${summary.max.toFixed(6)}
  极差: ${summary.range.toFixed(6)}
  方差: ${summary.variance.toFixed(6)}
  标准差: ${summary.stdDev.toFixed(6)}
  样本标准差: ${summary.sampleStdDev.toFixed(6)}
  变异系数: ${(summary.sampleStdDev/summary.mean*100).toFixed(2)}%

分位数:
  Q1 (25%): ${summary.q1.toFixed(6)}
  Q2 (50%): ${summary.q2.toFixed(6)}
  Q3 (75%): ${summary.q3.toFixed(6)}
  IQR: ${summary.iqr.toFixed(6)}

分布形态:
  偏度: ${summary.skewness.toFixed(6)}
  峰度: ${summary.kurtosis.toFixed(6)}`;
      
      setStatsResult(result);
      saveToHistory('统计分析', `n=${summary.count}, μ=${summary.mean.toFixed(4)}`, 'stats');
    } catch (e: any) {
      setStatsResult(`错误: ${e.message}`);
    }
  };

  // ==================== 设置 ====================
  
  const updateSetting = <K extends keyof Settings>(key: K, value: Settings[K]) => {
    const newSettings = { ...settings, [key]: value };
    setSettings(newSettings);
    storage.saveSettings(newSettings);
    if (key === 'angleUnit') calculator.setAngleUnit(value as Settings['angleUnit']);
  };

  // ==================== 渲染 ====================

  return (
    <div className={`app ${settings.theme}`}>
      <header className="header">
        <h1>测绘计算器</h1>
        <span className="version">Pro</span>
      </header>

      <main className="main-content">
        {/* 首页 */}
        {activeTab === 'home' && (
          <div className="home-page">
            <div className="welcome">
              <h2>专业测绘计算器</h2>
              <p>科学计算 · 测绘计算 · 统计分析</p>
            </div>
            <div className="quick-grid">
              <button onClick={() => setActiveTab('calc')}><span>🔢</span><label>科学计算</label></button>
              <button onClick={() => setActiveTab('survey')}><span>📐</span><label>测绘计算</label></button>
              <button onClick={() => setActiveTab('stats')}><span>📊</span><label>统计分析</label></button>
              <button onClick={() => setActiveTab('settings')}><span>⚙️</span><label>设置</label></button>
            </div>
            <div className="history-section">
              <h3>计算历史</h3>
              {history.slice(0, 10).map(item => (
                <div key={item.id} className="history-row" onClick={() => setDisplay(item.result)}>
                  <span className="expr">{item.expression}</span>
                  <span className="result">{item.result}</span>
                </div>
              ))}
              {history.length === 0 && <p className="empty">暂无历史记录</p>}
            </div>
          </div>
        )}

        {/* 科学计算器 */}
        {activeTab === 'calc' && (
          <div className="calc-page">
            <div className="calc-display">
              <div className="expr">{expression}</div>
              <div className="result">{display}</div>
              {hasMemory && <div className="memory-indicator">M</div>}
            </div>
            
            <div className="calc-toolbar">
              <div className="angle-unit">
                {(['DEG', 'RAD', 'GRAD'] as const).map(u => (
                  <button key={u} className={settings.angleUnit === u ? 'active' : ''} 
                    onClick={() => updateSetting('angleUnit', u)}>{u}</button>
                ))}
              </div>
              <button className={showScientific ? 'active' : ''} onClick={() => setShowScientific(!showScientific)}>
                {showScientific ? '简化' : '科学'}
              </button>
            </div>

            <div className="calc-buttons">
              {showScientific && (
                <div className="sci-panel">
                  <div className="sci-row">
                    <button onClick={() => insertFunction('sin')}>sin</button>
                    <button onClick={() => insertFunction('cos')}>cos</button>
                    <button onClick={() => insertFunction('tan')}>tan</button>
                    <button onClick={() => insertFunction('sinh')}>sinh</button>
                    <button onClick={() => insertFunction('cosh')}>cosh</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={() => insertFunction('asin')}>asin</button>
                    <button onClick={() => insertFunction('acos')}>acos</button>
                    <button onClick={() => insertFunction('atan')}>atan</button>
                    <button onClick={() => insertFunction('tanh')}>tanh</button>
                    <button onClick={() => applyFunction('π')}>π</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={() => insertFunction('ln')}>ln</button>
                    <button onClick={() => insertFunction('log')}>log</button>
                    <button onClick={() => insertFunction('log2')}>log₂</button>
                    <button onClick={() => applyFunction('e')}>e</button>
                    <button onClick={() => applyFunction('eˣ')}>eˣ</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={() => insertFunction('√')}>√</button>
                    <button onClick={() => insertFunction('∛')}>∛</button>
                    <button onClick={() => applyFunction('x²')}>x²</button>
                    <button onClick={() => applyFunction('x³')}>x³</button>
                    <button onClick={() => appendToDisplay('^')}>^</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={() => applyFunction('n!')}>n!</button>
                    <button onClick={() => applyFunction('1/x')}>1/x</button>
                    <button onClick={() => insertFunction('abs')}>|x|</button>
                    <button onClick={() => appendToDisplay('(')}>(</button>
                    <button onClick={() => appendToDisplay(')')}>)</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={memClear}>MC</button>
                    <button onClick={memRecall}>MR</button>
                    <button onClick={memAdd}>M+</button>
                    <button onClick={memSub}>M-</button>
                    <button onClick={() => applyFunction('rand')}>Rnd</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={() => insertDMS('°')}>°</button>
                    <button onClick={() => insertDMS("'")}>′</button>
                    <button onClick={() => insertDMS('"')}>″</button>
                    <button onClick={deg2dms}>D→DMS</button>
                    <button onClick={dms2deg}>DMS→D</button>
                  </div>
                  <div className="sci-row">
                    <button onClick={deg2rad}>D→R</button>
                    <button onClick={rad2deg}>R→D</button>
                    <button onClick={() => applyFunction('floor')}>⌊x⌋</button>
                    <button onClick={() => applyFunction('ceil')}>⌈x⌉</button>
                    <button onClick={() => applyFunction('ANS')}>ANS</button>
                  </div>
                </div>
              )}
              
              <div className="num-panel">
                <div className="num-row">
                  <button className="clear" onClick={clearAll}>AC</button>
                  <button className="func" onClick={toggleSign}>±</button>
                  <button className="func" onClick={() => applyFunction('%')}>%</button>
                  <button className="op" onClick={() => appendToDisplay('÷')}>÷</button>
                </div>
                <div className="num-row">
                  <button onClick={() => appendToDisplay('7')}>7</button>
                  <button onClick={() => appendToDisplay('8')}>8</button>
                  <button onClick={() => appendToDisplay('9')}>9</button>
                  <button className="op" onClick={() => appendToDisplay('×')}>×</button>
                </div>
                <div className="num-row">
                  <button onClick={() => appendToDisplay('4')}>4</button>
                  <button onClick={() => appendToDisplay('5')}>5</button>
                  <button onClick={() => appendToDisplay('6')}>6</button>
                  <button className="op" onClick={() => appendToDisplay('-')}>−</button>
                </div>
                <div className="num-row">
                  <button onClick={() => appendToDisplay('1')}>1</button>
                  <button onClick={() => appendToDisplay('2')}>2</button>
                  <button onClick={() => appendToDisplay('3')}>3</button>
                  <button className="op" onClick={() => appendToDisplay('+')}>+</button>
                </div>
                <div className="num-row">
                  <button onClick={() => appendToDisplay('0')}>0</button>
                  <button onClick={() => appendToDisplay('.')}>.</button>
                  <button className="func" onClick={backspace}>⌫</button>
                  <button className="equals" onClick={calculate}>=</button>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* 测绘计算 */}
        {activeTab === 'survey' && (
          <div className="survey-page">
            <div className="survey-types">
              {[
                { key: 'forward', label: '坐标正算' },
                { key: 'inverse', label: '坐标反算' },
                { key: 'forward_intersection', label: '前方交会' },
                { key: 'resection', label: '后方交会' },
                { key: 'side_shot', label: '侧方交会' },
                { key: 'area', label: '面积计算' },
                { key: 'gauss_forward', label: '高斯正算' },
                { key: 'gauss_inverse', label: '高斯反算' },
                { key: 'curve', label: '曲线要素' },
                { key: 'traverse', label: '导线平差' },
                { key: 'transform4', label: '四参数' },
                { key: 'leveling', label: '水准测量' },
                { key: 'earthwork', label: '土方计算' },
              ].map(({ key, label }) => (
                <button key={key} className={surveyType === key ? 'active' : ''} 
                  onClick={() => { setSurveyType(key); setSurveyResult(''); setSurveyInputs({}); }}>{label}</button>
              ))}
            </div>
            
            <div className="survey-form">
              {surveyType === 'forward' && (
                <>
                  <div className="input-group"><label>已知点X</label><input type="number" value={surveyInputs.x0||''} onChange={e=>handleSurveyInput('x0',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>已知点Y</label><input type="number" value={surveyInputs.y0||''} onChange={e=>handleSurveyInput('y0',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>方位角(°)</label><input type="number" value={surveyInputs.azimuth||''} onChange={e=>handleSurveyInput('azimuth',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>距离(m)</label><input type="number" value={surveyInputs.distance||''} onChange={e=>handleSurveyInput('distance',e.target.value)} placeholder="100"/></div>
                </>
              )}
              {surveyType === 'inverse' && (
                <>
                  <div className="input-group"><label>起点X</label><input type="number" value={surveyInputs.x1||''} onChange={e=>handleSurveyInput('x1',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>起点Y</label><input type="number" value={surveyInputs.y1||''} onChange={e=>handleSurveyInput('y1',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>终点X</label><input type="number" value={surveyInputs.x2||''} onChange={e=>handleSurveyInput('x2',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>终点Y</label><input type="number" value={surveyInputs.y2||''} onChange={e=>handleSurveyInput('y2',e.target.value)} placeholder="100"/></div>
                </>
              )}
              {surveyType === 'forward_intersection' && (
                <>
                  <div className="input-group"><label>A点X</label><input type="number" value={surveyInputs.xa||''} onChange={e=>handleSurveyInput('xa',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>A点Y</label><input type="number" value={surveyInputs.ya||''} onChange={e=>handleSurveyInput('ya',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>B点X</label><input type="number" value={surveyInputs.xb||''} onChange={e=>handleSurveyInput('xb',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>B点Y</label><input type="number" value={surveyInputs.yb||''} onChange={e=>handleSurveyInput('yb',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>∠PAB(°)</label><input type="number" value={surveyInputs.angleA||''} onChange={e=>handleSurveyInput('angleA',e.target.value)} placeholder="45"/></div>
                  <div className="input-group"><label>∠PBA(°)</label><input type="number" value={surveyInputs.angleB||''} onChange={e=>handleSurveyInput('angleB',e.target.value)} placeholder="45"/></div>
                </>
              )}
              {surveyType === 'resection' && (
                <>
                  <div className="input-group"><label>A点X</label><input type="number" value={surveyInputs.xa||''} onChange={e=>handleSurveyInput('xa',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>A点Y</label><input type="number" value={surveyInputs.ya||''} onChange={e=>handleSurveyInput('ya',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>B点X</label><input type="number" value={surveyInputs.xb||''} onChange={e=>handleSurveyInput('xb',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>B点Y</label><input type="number" value={surveyInputs.yb||''} onChange={e=>handleSurveyInput('yb',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>C点X</label><input type="number" value={surveyInputs.xc||''} onChange={e=>handleSurveyInput('xc',e.target.value)} placeholder="50"/></div>
                  <div className="input-group"><label>C点Y</label><input type="number" value={surveyInputs.yc||''} onChange={e=>handleSurveyInput('yc',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>∠APB(°)</label><input type="number" value={surveyInputs.alpha||''} onChange={e=>handleSurveyInput('alpha',e.target.value)} placeholder="60"/></div>
                  <div className="input-group"><label>∠BPC(°)</label><input type="number" value={surveyInputs.beta||''} onChange={e=>handleSurveyInput('beta',e.target.value)} placeholder="60"/></div>
                </>
              )}
              {surveyType === 'side_shot' && (
                <>
                  <div className="input-group"><label>测站X</label><input type="number" value={surveyInputs.x0||''} onChange={e=>handleSurveyInput('x0',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>测站Y</label><input type="number" value={surveyInputs.y0||''} onChange={e=>handleSurveyInput('y0',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>后视方位角(°)</label><input type="number" value={surveyInputs.backAzimuth||''} onChange={e=>handleSurveyInput('backAzimuth',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>水平角(°)</label><input type="number" value={surveyInputs.angle||''} onChange={e=>handleSurveyInput('angle',e.target.value)} placeholder="90"/></div>
                  <div className="input-group"><label>距离(m)</label><input type="number" value={surveyInputs.distance||''} onChange={e=>handleSurveyInput('distance',e.target.value)} placeholder="100"/></div>
                </>
              )}
              {surveyType === 'area' && (
                <div className="input-group full">
                  <label>多边形顶点(X,Y每行一个)</label>
                  <textarea value={surveyInputs.points||'0,0\n100,0\n100,100\n0,100'} onChange={e=>handleSurveyInput('points',e.target.value)} rows={6}/>
                </div>
              )}
              {surveyType === 'gauss_forward' && (
                <>
                  <div className="input-group"><label>纬度(°)</label><input type="number" value={surveyInputs.lat||''} onChange={e=>handleSurveyInput('lat',e.target.value)} placeholder="30"/></div>
                  <div className="input-group"><label>经度(°)</label><input type="number" value={surveyInputs.lon||''} onChange={e=>handleSurveyInput('lon',e.target.value)} placeholder="120"/></div>
                  <div className="input-group"><label>中央子午线(可选)</label><input type="number" value={surveyInputs.cm||''} onChange={e=>handleSurveyInput('cm',e.target.value)} placeholder="自动计算"/></div>
                </>
              )}
              {surveyType === 'gauss_inverse' && (
                <>
                  <div className="input-group"><label>X坐标(m)</label><input type="number" value={surveyInputs.gx||''} onChange={e=>handleSurveyInput('gx',e.target.value)} placeholder="3000000"/></div>
                  <div className="input-group"><label>Y坐标(m)</label><input type="number" value={surveyInputs.gy||''} onChange={e=>handleSurveyInput('gy',e.target.value)} placeholder="500000"/></div>
                  <div className="input-group"><label>中央子午线(°)</label><input type="number" value={surveyInputs.cm||''} onChange={e=>handleSurveyInput('cm',e.target.value)} placeholder="120"/></div>
                </>
              )}
              {surveyType === 'curve' && (
                <>
                  <div className="input-group"><label>半径R(m)</label><input type="number" value={surveyInputs.radius||''} onChange={e=>handleSurveyInput('radius',e.target.value)} placeholder="500"/></div>
                  <div className="input-group"><label>转向角(°)</label><input type="number" value={surveyInputs.deflection||''} onChange={e=>handleSurveyInput('deflection',e.target.value)} placeholder="30"/></div>
                </>
              )}
              {surveyType === 'traverse' && (
                <>
                  <div className="input-group"><label>起点X</label><input type="number" value={surveyInputs.startX||''} onChange={e=>handleSurveyInput('startX',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>起点Y</label><input type="number" value={surveyInputs.startY||''} onChange={e=>handleSurveyInput('startY',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>起始方位角(°)</label><input type="number" value={surveyInputs.startAz||''} onChange={e=>handleSurveyInput('startAz',e.target.value)} placeholder="0"/></div>
                  <div className="input-group full"><label>测站数据(角度,距离 每行一站)</label><textarea value={surveyInputs.stations||'90,100\n90,100\n90,100\n90,100'} onChange={e=>handleSurveyInput('stations',e.target.value)} rows={5}/></div>
                </>
              )}
              {surveyType === 'transform4' && (
                <>
                  <div className="input-group"><label>平移ΔX(m)</label><input type="number" value={surveyInputs.dx||''} onChange={e=>handleSurveyInput('dx',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>平移ΔY(m)</label><input type="number" value={surveyInputs.dy||''} onChange={e=>handleSurveyInput('dy',e.target.value)} placeholder="200"/></div>
                  <div className="input-group"><label>尺度因子</label><input type="number" value={surveyInputs.scale||''} onChange={e=>handleSurveyInput('scale',e.target.value)} placeholder="1" step="0.0000001"/></div>
                  <div className="input-group"><label>旋转角(°)</label><input type="number" value={surveyInputs.rotation||''} onChange={e=>handleSurveyInput('rotation',e.target.value)} placeholder="0"/></div>
                  <div className="input-group"><label>原X坐标</label><input type="number" value={surveyInputs.tx||''} onChange={e=>handleSurveyInput('tx',e.target.value)} placeholder="1000"/></div>
                  <div className="input-group"><label>原Y坐标</label><input type="number" value={surveyInputs.ty||''} onChange={e=>handleSurveyInput('ty',e.target.value)} placeholder="2000"/></div>
                </>
              )}
              {surveyType === 'leveling' && (
                <div className="input-group full">
                  <label>水准数据(第一行起始高程，后续每行:后视,前视)</label>
                  <textarea value={surveyInputs.heights||'100.000\n1.234,-2.345\n1.567,-1.890\n1.123,-2.456'} onChange={e=>handleSurveyInput('heights',e.target.value)} rows={6}/>
                </div>
              )}
              {surveyType === 'earthwork' && (
                <>
                  <div className="input-group"><label>断面1面积(m²)</label><input type="number" value={surveyInputs.area1||''} onChange={e=>handleSurveyInput('area1',e.target.value)} placeholder="100"/></div>
                  <div className="input-group"><label>断面2面积(m²)</label><input type="number" value={surveyInputs.area2||''} onChange={e=>handleSurveyInput('area2',e.target.value)} placeholder="120"/></div>
                  <div className="input-group"><label>断面间距(m)</label><input type="number" value={surveyInputs.edist||''} onChange={e=>handleSurveyInput('edist',e.target.value)} placeholder="50"/></div>
                </>
              )}
              
              <button className="calc-btn" onClick={calculateSurvey}>计 算</button>
            </div>
            
            {surveyResult && (
              <div className="survey-result">
                <pre>{surveyResult}</pre>
              </div>
            )}
          </div>
        )}

        {/* 统计分析 */}
        {activeTab === 'stats' && (
          <div className="stats-page">
            <div className="stats-input">
              <label>输入数据(逗号、空格或换行分隔)</label>
              <textarea value={statsData} onChange={e => setStatsData(e.target.value)} 
                placeholder="1, 2, 3, 4, 5&#10;或每行一个数据" rows={6}/>
              <button className="calc-btn" onClick={calculateStats}>统计分析</button>
            </div>
            {statsResult && (
              <div className="stats-result">
                <pre>{statsResult}</pre>
              </div>
            )}
          </div>
        )}

        {/* 设置 */}
        {activeTab === 'settings' && (
          <div className="settings-page">
            <div className="settings-group">
              <h3>显示设置</h3>
              <div className="setting-item">
                <span>主题</span>
                <div className="setting-options">
                  <button className={settings.theme==='dark'?'active':''} onClick={()=>updateSetting('theme','dark')}>深色</button>
                  <button className={settings.theme==='light'?'active':''} onClick={()=>updateSetting('theme','light')}>浅色</button>
                </div>
              </div>
              <div className="setting-item">
                <span>角度单位</span>
                <div className="setting-options">
                  {(['DEG','RAD','GRAD'] as const).map(u=>(
                    <button key={u} className={settings.angleUnit===u?'active':''} onClick={()=>updateSetting('angleUnit',u)}>{u}</button>
                  ))}
                </div>
              </div>
              <div className="setting-item">
                <span>小数位数</span>
                <div className="setting-options">
                  <input type="range" min="1" max="15" value={settings.precision} onChange={e=>updateSetting('precision',parseInt(e.target.value))}/>
                  <span>{settings.precision}</span>
                </div>
              </div>
              <div className="setting-item">
                <span>震动反馈</span>
                <button className={`toggle ${settings.vibration?'active':''}`} onClick={()=>updateSetting('vibration',!settings.vibration)}/>
              </div>
            </div>
            <div className="settings-group">
              <h3>数据管理</h3>
              <button className="action-btn" onClick={()=>{storage.clearHistory();setHistory([]);}}>清除历史记录</button>
              <button className="action-btn" onClick={()=>alert(storage.exportAll())}>导出数据</button>
            </div>
            <div className="app-info">
              <div className="logo">📐</div>
              <div className="name">测绘计算器 Pro</div>
              <div className="ver">版本 2.0.0</div>
            </div>
          </div>
        )}
      </main>

      <nav className="bottom-nav">
        {[
          { key: 'home', icon: '🏠', label: '首页' },
          { key: 'calc', icon: '🔢', label: '计算' },
          { key: 'survey', icon: '📐', label: '测绘' },
          { key: 'stats', icon: '📊', label: '统计' },
          { key: 'settings', icon: '⚙️', label: '设置' },
        ].map(({ key, icon, label }) => (
          <button key={key} className={`nav-item ${activeTab === key ? 'active' : ''}`}
            onClick={() => setActiveTab(key as TabType)}>
            <span className="icon">{icon}</span>
            <span className="label">{label}</span>
          </button>
        ))}
      </nav>
    </div>
  );
}

export default App;
