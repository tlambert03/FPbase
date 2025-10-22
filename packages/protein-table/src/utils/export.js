/**
 * Convert array of objects to CSV string
 */
function arrayToCSV(data) {
  if (data.length === 0) return ""

  // Get headers from first object
  const headers = Object.keys(data[0])

  // Create CSV rows
  const rows = data.map((row) => {
    return headers.map((header) => {
      const value = row[header]
      // Escape quotes and wrap in quotes if contains comma, quote, or newline
      if (value == null) return ""
      const stringValue = String(value)
      if (stringValue.includes(",") || stringValue.includes('"') || stringValue.includes("\n")) {
        return `"${stringValue.replace(/"/g, '""')}"`
      }
      return stringValue
    }).join(",")
  })

  // Combine headers and rows
  return [headers.join(","), ...rows].join("\n")
}

/**
 * Export table data to CSV format
 */
export function exportToCSV(data, filename = "proteins.csv") {
  const csv = arrayToCSV(data)
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" })

  // Create download link
  const link = document.createElement("a")
  const url = URL.createObjectURL(blob)
  link.setAttribute("href", url)
  link.setAttribute("download", filename)
  link.style.visibility = "hidden"

  // Trigger download
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}

/**
 * Prepare row data for export (flatten states into separate rows)
 */
export function prepareExportData(proteins) {
  const rows = []

  proteins.forEach((protein) => {
    const states = protein.states?.length > 0 ? protein.states : [{}]

    states.forEach((state) => {
      // Skip dark states
      if (state.is_dark) return

      const row = {
        Name: protein.name,
        "State": state.name !== "default" ? state.name : "",
        "Ex max (nm)": state.ex_max || "",
        "Em max (nm)": state.em_max || "",
        "Stokes Shift (nm)": state.stokes || "",
        "Extinction Coefficient": state.ext_coeff || "",
        "Quantum Yield": state.qy || "",
        "Brightness": state.brightness || "",
        "pKa": state.pka || "",
        "Oligomerization": protein.agg || "",
        "Maturation (min)": state.maturation || "",
        "Lifetime (ns)": state.lifetime || "",
        "Molecular Weight (kDa)": protein.weight || "",
        "Year": protein.year || "",
        "Switch Type": protein.switch_type || "",
        "Aliases": Array.isArray(protein.aliases) ? protein.aliases.join(", ") : "",
      }
      rows.push(row)
    })
  })

  return rows
}
