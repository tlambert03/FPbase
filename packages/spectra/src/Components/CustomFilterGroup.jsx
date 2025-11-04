import AddIcon from "@mui/icons-material/Add"
import DeleteIcon from "@mui/icons-material/Delete"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import IconButton from "@mui/material/IconButton"
import { useEffect, useRef, useState } from "react"
import { useSpectraStore } from "../store/spectraStore"
import CustomFilterCreator from "./CustomFilterCreator"
import { categoryIcon } from "./FaIcon"

const CustomFilterGroup = ({ activeSpectra, showAddButton = true }) => {
  const filterCounter = useRef(0)
  const [customFilters, setFilters] = useState([])
  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)
  const addCustomFilter = useSpectraStore((state) => state.addCustomFilter)
  const removeCustomFilter = useSpectraStore((state) => state.removeCustomFilter)

  const addRow = () => {
    const newId = `$cf${filterCounter.current++}`
    // Add to store with default params
    addCustomFilter(newId, { type: "BP", center: 525, width: 50, transmission: 90 })
    // Add to active spectra
    updateActiveSpectra([newId])
    // Add to local list
    setFilters([...customFilters, newId])
  }

  const removeRow = (filterId) => {
    // Remove from local list
    setFilters(customFilters.filter((id) => id !== filterId))
    // Remove from store
    removeCustomFilter(filterId)
    // Remove from active spectra
    updateActiveSpectra([], [filterId])
  }

  // Sync with activeSpectra (bidirectional: detect additions AND removals)
  useEffect(() => {
    const activeFilters = activeSpectra.filter((id) => id.startsWith("$cf"))

    // Detect new filters added from URL
    const newFilters = activeFilters.filter((id) => !customFilters.includes(id))
    if (newFilters.length) {
      const inds = newFilters.map((id) => Number.parseInt(id.replace("$cf", ""), 10))
      filterCounter.current = Math.max(...inds, filterCounter.current) + 1
    }

    // Detect removed filters (e.g., from clearAllSpectra)
    const removedFilters = customFilters.filter((id) => !activeFilters.includes(id))

    // Update local state if there are changes
    if (newFilters.length > 0 || removedFilters.length > 0) {
      setFilters(activeFilters)
    }
  }, [activeSpectra, customFilters]) // eslint-disable-line

  return (
    <div>
      {customFilters.sort().map((filter) => (
        <div style={{ width: "100%", margin: "4px 0" }} key={filter}>
          <Box display="flex" alignItems="center">
            {categoryIcon("CF", "rgba(0,0,50,0.4)", {
              style: {
                position: "relative",
                top: 0,
                left: 2,
                height: "1.3rem",
                marginRight: 10,
              },
            })}
            <Box flexGrow={1}>
              <CustomFilterCreator key={filter} id={filter} />
            </Box>
            <Box>
              <IconButton
                aria-label="Delete"
                color="secondary"
                tabIndex={-1}
                onClick={() => removeRow(filter)}
                style={{
                  padding: "6px 6px",
                  marginRight: 2,
                  marginLeft: 2,
                }}
              >
                <DeleteIcon />
              </IconButton>
            </Box>
          </Box>
        </div>
      ))}
      {showAddButton && (
        <Button
          variant="contained"
          color="primary"
          onClick={() => addRow()}
          style={{ marginTop: 8, marginLeft: 34 }}
        >
          <AddIcon />
          {`Add Custom Filter`}
        </Button>
      )}
    </div>
  )
}

export default CustomFilterGroup
