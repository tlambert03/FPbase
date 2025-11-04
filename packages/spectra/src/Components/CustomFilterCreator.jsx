import Box from "@mui/material/Box"
import ToggleButton from "@mui/material/ToggleButton"
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup"
import Tooltip from "@mui/material/Tooltip"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import React from "react"
import { useSpectraStore } from "../store/spectraStore"
import InputSlider from "./InputSlider"

export const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444",
  },
})

const CustomFilterCreator = React.memo(function CustomFilterCreator({ id }) {
  const classes = useStyles()

  // Read params directly from store (single source of truth)
  const params = useSpectraStore((state) => state.customFilters[id])
  const updateCustomFilter = useSpectraStore((state) => state.updateCustomFilter)

  // Default values if params don't exist yet
  const type = params?.type || "BP"
  const center = params?.center || 525
  const width = params?.width || 50
  const transmission = params?.transmission || 90

  const handleType = (_event, newType) => {
    if (newType) {
      updateCustomFilter(id, { type: newType })
    }
  }

  return (
    <div
      style={{
        padding: "10px 10px 0px 10px",
        border: "1px solid #ccc",
        borderRadius: 4,
      }}
    >
      <Box display="flex" flexWrap="wrap">
        <Typography style={{ margin: "8px 10px 3px" }}>Custom Filter</Typography>
        <ToggleButtonGroup
          size="small"
          value={type}
          exclusive
          onChange={handleType}
          style={{ marginBottom: 10, marginLeft: 10 }}
        >
          <Tooltip title="Longpass Filter">
            <ToggleButton value="LP">LP</ToggleButton>
          </Tooltip>
          <Tooltip title="Shortpass Filter">
            <ToggleButton value="SP">SP</ToggleButton>
          </Tooltip>
          <Tooltip title="Bandpass Filter">
            <ToggleButton value="BP">BP</ToggleButton>
          </Tooltip>
        </ToggleButtonGroup>

        <Box flexGrow={2}>
          <div style={{ margin: "0 12px", minWidth: 180 }}>
            <Typography className={classes.label}>
              {type === "BP" ? "Center Wavelength" : "Edge"}
            </Typography>
            <InputSlider
              value={center}
              setValue={(value) => updateCustomFilter(id, { center: value })}
              valueLabelDisplay="auto"
              aria-labelledby="range-slider"
            />
          </div>
        </Box>
        {type === "BP" && (
          <Box flexGrow={1}>
            <div
              style={{
                margin: "0 10px",
                minWidth: 150,
              }}
            >
              <Typography className={classes.label}>Bandwidth</Typography>
              <InputSlider
                value={width}
                setValue={(value) => updateCustomFilter(id, { width: value })}
                min={1}
                max={300}
                valueLabelDisplay="auto"
                aria-labelledby="range-slider"
              />
            </div>
          </Box>
        )}
        <Box flexGrow={1}>
          <div
            style={{
              margin: "0 10px",
              minWidth: 160,
            }}
          >
            <Typography className={classes.label}>Transmission %</Typography>
            <InputSlider
              value={transmission}
              setValue={(value) => updateCustomFilter(id, { transmission: value })}
              min={1}
              max={100}
              valueLabelDisplay="auto"
              aria-labelledby="range-slider"
            />
          </div>
        </Box>
      </Box>
    </div>
  )
})

export default CustomFilterCreator
