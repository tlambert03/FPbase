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
  plugins: [react()],
})
