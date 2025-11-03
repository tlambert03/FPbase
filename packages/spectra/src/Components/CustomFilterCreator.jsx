import Box from "@mui/material/Box"
import ToggleButton from "@mui/material/ToggleButton"
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import React, { useEffect } from "react"
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

  const [filterID, _type, _center, _width, _trans] = id.split("_")
  const [type, setType] = React.useState((_type || "").toUpperCase() || "BP")
  const [center, setCenter] = React.useState(_center || 525)
  const [width, setWidth] = React.useState(_width || 50)
  const [trans, setTrans] = React.useState(_trans || 90)

  // Track the previous full ID to properly remove it when parameters change
  const prevIdRef = React.useRef(id)

  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)
  useEffect(() => {
    const newId = `${filterID}_${type}_${center}_${width}_${trans}`
    const oldId = prevIdRef.current

    // Only update if the ID actually changed
    if (newId !== oldId) {
      updateActiveSpectra([newId], [oldId])
      prevIdRef.current = newId
    }
  }, [width, center, type, trans, filterID, updateActiveSpectra])

  const handleType = (_event, newType) => {
    if (newType) {
      setType(newType)
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
        <Typography style={{ margin: "8px 10px 3px" }}>Custom Filter Type</Typography>
        <ToggleButtonGroup
          size="small"
          value={type}
          exclusive
          onChange={handleType}
          style={{ marginBottom: 10, marginLeft: 10 }}
        >
          <ToggleButton value="LP">LP</ToggleButton>
          <ToggleButton value="SP">SP</ToggleButton>
          <ToggleButton value="BP">BP</ToggleButton>
        </ToggleButtonGroup>

        <Box flexGrow={2}>
          <div style={{ margin: "0 12px", minWidth: 180 }}>
            <Typography className={classes.label}>
              {type === "BP" ? "Center Wavelength" : "Edge"}
            </Typography>
            <InputSlider
              value={center}
              setValue={setCenter}
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
                setValue={setWidth}
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
              value={trans}
              setValue={setTrans}
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
