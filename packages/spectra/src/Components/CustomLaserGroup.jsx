import AddIcon from "@mui/icons-material/Add"
import DeleteIcon from "@mui/icons-material/Delete"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import IconButton from "@mui/material/IconButton"
import React, { useEffect, useRef, useState } from "react"
import { useSpectraStore } from "../store/spectraStore"
import CustomLaserCreator from "./CustomLaserCreator"
import { categoryIcon } from "./FaIcon"

const CustomLaserGroup = React.memo(function CustomLaserGroup({
  activeSpectra,
  showAddButton = true,
}) {
  const laserCounter = useRef(0)
  const [customLasers, setLasers] = useState([])
  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)
  const addCustomLaser = useSpectraStore((state) => state.addCustomLaser)
  const removeCustomLaser = useSpectraStore((state) => state.removeCustomLaser)
  const exNorm = useSpectraStore((state) => state.exNorm)
  const setExNormStore = useSpectraStore((state) => state.setExNorm)
  const normID = exNorm?.[1]

  const setExNorm = React.useCallback((data) => setExNormStore(data), [setExNormStore])

  const clearNorm = React.useCallback(() => setExNormStore([null, null]), [setExNormStore])

  // Sync with activeSpectra to detect lasers added from URL
  useEffect(() => {
    if (activeSpectra && activeSpectra.length > 0) {
      const newLasers = activeSpectra.filter(
        (id) => id.startsWith("$cl") && !customLasers.includes(id)
      )
      if (newLasers.length) {
        const inds = newLasers.map((id) => Number.parseInt(id.replace("$cl", ""), 10))
        laserCounter.current = Math.max(...inds, laserCounter.current) + 1
        setLasers([...customLasers, ...newLasers])
      }
    }
  }, [activeSpectra, customLasers]) // eslint-disable-line

  const addRow = () => {
    const newId = `$cl${laserCounter.current++}`
    // Add to store with default params
    addCustomLaser(newId, { wavelength: 488 })
    // Add to active spectra
    updateActiveSpectra([newId])
    // Add to local list
    setLasers([...customLasers, newId])
  }

  const removeRow = (laserId) => {
    if (laserId === normID) {
      clearNorm()
    }
    // Remove from local list
    setLasers(customLasers.filter((id) => id !== laserId))
    // Remove from store
    removeCustomLaser(laserId)
    // Remove from active spectra
    updateActiveSpectra([], [laserId])
  }

  return (
    <div>
      {customLasers.sort().map((laser) => (
        <div style={{ width: "100%", margin: "4px 0" }} key={laser}>
          <Box display="flex" alignItems="center">
            {categoryIcon("CL", "rgba(0,0,50,0.4)", {
              style: {
                position: "relative",
                top: 0,
                left: 2,
                height: "1.3rem",
                marginRight: 10,
              },
            })}
            <Box flexGrow={1}>
              <CustomLaserCreator
                key={laser}
                id={laser}
                setExNorm={setExNorm}
                clearNorm={clearNorm}
                normID={normID}
              />
            </Box>
            <Box>
              <IconButton
                aria-label="Delete"
                color="secondary"
                tabIndex={-1}
                onClick={() => removeRow(laser)}
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
          {`Add Laser`}
        </Button>
      )}
    </div>
  )
})

export default CustomLaserGroup
