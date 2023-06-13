import { defineConfig } from "vite"
import react from "@vitejs/plugin-react";

import pkg from './package.json';

// const external = Object.keys(pkg.dependencies);
// external.push('react', 'react-dom')
const external = ['react']
export default defineConfig({
  build: {
    lib: {
      entry: 'src/index.jsx',
      formats: ['es'],
    },
    rollupOptions: { external }
  },
  server: {
    proxy: {
      '/graphql': 'http://127.0.0.1:8000',
      '/api': 'http://127.0.0.1:8000',
    },
  },
  plugins: [react()],
})
