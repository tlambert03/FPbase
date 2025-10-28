import { StrictMode } from "react"
import { createRoot } from "react-dom/client"
import App from "./App"

export default function mount(el) {
  const root = createRoot(el)
  root.render(
    <StrictMode>
      <App />
    </StrictMode>
  )
}
