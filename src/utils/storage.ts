/**
 * 本地存储管理
 */

export interface HistoryItem {
  id: string;
  expression: string;
  result: string;
  timestamp: number;
  type: 'calc' | 'survey' | 'stats' | 'matrix';
}

export interface Project {
  id: string;
  name: string;
  type: string;
  data: any;
  createdAt: number;
  updatedAt: number;
}

export interface Settings {
  theme: 'light' | 'dark';
  angleUnit: 'DEG' | 'RAD' | 'GRAD';
  precision: number;
  vibration: boolean;
  soundEnabled: boolean;
  scientificNotation: boolean;
}

const STORAGE_KEYS = {
  HISTORY: 'survey_calc_history',
  PROJECTS: 'survey_calc_projects',
  SETTINGS: 'survey_calc_settings',
  VARIABLES: 'survey_calc_variables',
  MATRICES: 'survey_calc_matrices',
};

class StorageManager {
  // 设置
  getSettings(): Settings {
    const data = localStorage.getItem(STORAGE_KEYS.SETTINGS);
    return data ? JSON.parse(data) : {
      theme: 'dark',
      angleUnit: 'DEG',
      precision: 8,
      vibration: true,
      soundEnabled: false,
      scientificNotation: false,
    };
  }

  saveSettings(settings: Settings): void {
    localStorage.setItem(STORAGE_KEYS.SETTINGS, JSON.stringify(settings));
  }

  // 历史记录
  getHistory(limit: number = 100): HistoryItem[] {
    const data = localStorage.getItem(STORAGE_KEYS.HISTORY);
    const history: HistoryItem[] = data ? JSON.parse(data) : [];
    return history.slice(0, limit);
  }

  addHistory(item: Omit<HistoryItem, 'id' | 'timestamp'>): void {
    const history = this.getHistory(500);
    history.unshift({
      ...item,
      id: Date.now().toString(36) + Math.random().toString(36).substr(2, 9),
      timestamp: Date.now(),
    });
    localStorage.setItem(STORAGE_KEYS.HISTORY, JSON.stringify(history.slice(0, 500)));
  }

  clearHistory(): void {
    localStorage.setItem(STORAGE_KEYS.HISTORY, JSON.stringify([]));
  }

  deleteHistoryItem(id: string): void {
    const history = this.getHistory(500).filter(h => h.id !== id);
    localStorage.setItem(STORAGE_KEYS.HISTORY, JSON.stringify(history));
  }

  // 项目
  getProjects(): Project[] {
    const data = localStorage.getItem(STORAGE_KEYS.PROJECTS);
    return data ? JSON.parse(data) : [];
  }

  saveProject(project: Omit<Project, 'id' | 'createdAt' | 'updatedAt'>): string {
    const projects = this.getProjects();
    const id = Date.now().toString(36) + Math.random().toString(36).substr(2, 9);
    const now = Date.now();
    projects.unshift({
      ...project,
      id,
      createdAt: now,
      updatedAt: now,
    });
    localStorage.setItem(STORAGE_KEYS.PROJECTS, JSON.stringify(projects));
    return id;
  }

  updateProject(id: string, data: Partial<Project>): void {
    const projects = this.getProjects().map(p => 
      p.id === id ? { ...p, ...data, updatedAt: Date.now() } : p
    );
    localStorage.setItem(STORAGE_KEYS.PROJECTS, JSON.stringify(projects));
  }

  deleteProject(id: string): void {
    const projects = this.getProjects().filter(p => p.id !== id);
    localStorage.setItem(STORAGE_KEYS.PROJECTS, JSON.stringify(projects));
  }

  // 变量
  getVariables(): Map<string, number> {
    const data = localStorage.getItem(STORAGE_KEYS.VARIABLES);
    return data ? new Map(JSON.parse(data)) : new Map();
  }

  saveVariables(vars: Map<string, number>): void {
    localStorage.setItem(STORAGE_KEYS.VARIABLES, JSON.stringify(Array.from(vars)));
  }

  // 矩阵
  getMatrices(): Map<string, number[][]> {
    const data = localStorage.getItem(STORAGE_KEYS.MATRICES);
    return data ? new Map(JSON.parse(data)) : new Map();
  }

  saveMatrices(matrices: Map<string, number[][]>): void {
    localStorage.setItem(STORAGE_KEYS.MATRICES, JSON.stringify(Array.from(matrices)));
  }

  // 导出所有数据
  exportAll(): string {
    return JSON.stringify({
      settings: this.getSettings(),
      history: this.getHistory(10000),
      projects: this.getProjects(),
      variables: Array.from(this.getVariables()),
      matrices: Array.from(this.getMatrices()),
      exportTime: new Date().toISOString(),
    }, null, 2);
  }

  // 导入数据
  importAll(jsonString: string): void {
    const data = JSON.parse(jsonString);
    if (data.settings) localStorage.setItem(STORAGE_KEYS.SETTINGS, JSON.stringify(data.settings));
    if (data.history) localStorage.setItem(STORAGE_KEYS.HISTORY, JSON.stringify(data.history));
    if (data.projects) localStorage.setItem(STORAGE_KEYS.PROJECTS, JSON.stringify(data.projects));
    if (data.variables) localStorage.setItem(STORAGE_KEYS.VARIABLES, JSON.stringify(data.variables));
    if (data.matrices) localStorage.setItem(STORAGE_KEYS.MATRICES, JSON.stringify(data.matrices));
  }

  // 清除所有数据
  clearAll(): void {
    Object.values(STORAGE_KEYS).forEach(key => localStorage.removeItem(key));
  }
}

export const storage = new StorageManager();
export default StorageManager;
