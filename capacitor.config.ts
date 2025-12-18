import { CapacitorConfig } from '@capacitor/cli';

const config: CapacitorConfig = {
  appId: 'com.survey.calculator',
  appName: '测绘计算器Pro',
  webDir: 'build',
  server: {
    androidScheme: 'https'
  }
};

export default config;
