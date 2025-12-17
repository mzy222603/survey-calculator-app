# 测绘计算器 - 获取APK指南

## 项目概述

这是一个功能完善的测绘计算器Android应用，包含：

- **科学计算器**：三角函数、对数、指数、阶乘、排列组合等
- **测绘计算**：坐标正反算、交会计算、导线平差、高斯投影、曲线要素等
- **统计分析**：描述统计、回归分析（线性/指数/对数/幂）
- **矩阵运算**：行列式、逆矩阵、特征值、LU分解等
- **完全离线运行**，无需网络连接

---

## 获取APK的方法

### 方法一：使用GitHub Actions自动构建（推荐）

1. **安装Git**
   - 下载：https://git-scm.com/download/win
   - 安装后重启电脑

2. **创建GitHub账号**
   - 访问：https://github.com
   - 注册一个免费账号

3. **上传项目到GitHub**
   ```bash
   cd survey-calc-app
   git init
   git add .
   git commit -m "Initial commit"
   git branch -M main
   git remote add origin https://github.com/你的用户名/survey-calculator.git
   git push -u origin main
   ```

4. **等待自动构建**
   - 推送后，GitHub Actions会自动开始构建
   - 约5-10分钟后，在 Actions 页面可以下载APK
   - 或者在 Releases 页面下载

5. **下载APK**
   - 进入 GitHub 仓库 → Actions → 最新的workflow → Artifacts
   - 下载 `survey-calculator-debug.zip`
   - 解压后得到 `app-debug.apk`

---

### 方法二：在线构建服务

#### 使用 Codemagic（免费）

1. 访问 https://codemagic.io
2. 用GitHub账号登录
3. 添加你的仓库
4. 选择 "React Native / Capacitor" 项目类型
5. 点击 "Start new build"
6. 构建完成后下载APK

#### 使用 Expo EAS Build

1. 安装 Expo CLI: `npm install -g expo-cli eas-cli`
2. 配置并构建

---

### 方法三：本地构建（需要Android SDK）

1. **安装Android Studio**
   - 下载：https://developer.android.com/studio
   - 安装时勾选 Android SDK

2. **设置环境变量**
   ```
   ANDROID_HOME = C:\Users\你的用户名\AppData\Local\Android\Sdk
   ```

3. **构建APK**
   ```bash
   cd survey-calc-app
   npm run build
   npx cap sync android
   cd android
   ./gradlew assembleDebug
   ```

4. **获取APK**
   - 位置：`android/app/build/outputs/apk/debug/app-debug.apk`

---

## 快速本地预览

在生成APK之前，你可以在浏览器中预览应用：

```bash
cd survey-calc-app
npm start
```

然后在浏览器中打开 http://localhost:3000 查看效果。

---

## 文件结构

```
survey-calc-app/
├── src/
│   ├── App.tsx          # 主应用组件
│   ├── App.css          # 样式文件
│   └── utils/
│       ├── calculator.ts    # 科学计算引擎
│       ├── survey.ts        # 测绘计算引擎
│       ├── statistics.ts    # 统计分析引擎
│       ├── matrix.ts        # 矩阵计算引擎
│       └── storage.ts       # 数据存储管理
├── android/             # Android原生项目
├── .github/
│   └── workflows/
│       └── build-android.yml  # GitHub Actions配置
└── package.json
```

---

## 常见问题

**Q: 为什么选择React + Capacitor而不是其他方案？**
A: 因为Python/Kivy不支持Python 3.14，Flutter未安装。React + Capacitor是最可行的方案，且支持GitHub Actions云端构建。

**Q: APK安装失败怎么办？**
A: 请在手机设置中允许"安装未知来源应用"。

**Q: 应用是否需要网络？**
A: 不需要，完全离线运行。

---

## 技术栈

- React 19 + TypeScript
- Ionic React UI组件
- Capacitor 8 跨平台框架
- mathjs 数学库

版本：1.0.0
