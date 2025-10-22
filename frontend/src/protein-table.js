import mount from "@fpbase/protein-table"

// Mount the protein table app when the DOM is ready
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", () => {
    const container = document.getElementById("protein-table")
    if (container) {
      mount(container)
    }
  })
} else {
  const container = document.getElementById("protein-table")
  if (container) {
    mount(container)
  }
}
