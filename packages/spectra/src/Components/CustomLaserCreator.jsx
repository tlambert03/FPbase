import RadioButtonChecked from "@mui/icons-material/RadioButtonChecked"
import RadioButtonUnchecked from "@mui/icons-material/RadioButtonUnchecked"
import Box from "@mui/material/Box"
import Checkbox from "@mui/material/Checkbox"
import FormControlLabel from "@mui/material/FormControlLabel"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import { memo, useEffect, useRef, useState } from "react"
import { useSpectraStore } from "../store/spectraStore"
import InputSlider from "./InputSlider"

const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444",
  },
})

// $cl1
const CustomLaserCreator = memo(function CustomLaserCreator({ id, normID, setExNorm, clearNorm }) {
  const classes = useStyles()

  const [laserID, _wave] = id.split("_")
  const [wave, setWave] = useState(_wave || 488)

  // Track the previous full ID to properly remove it when wavelength changes
  const prevIdRef = useRef(id)

  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)
  useEffect(() => {
    const newId = `${laserID}_${wave}`
    const oldId = prevIdRef.current

    // Only update if the ID actually changed
    if (newId !== oldId) {
      updateActiveSpectra([newId], [oldId])
      prevIdRef.current = newId
    }

    if (laserID === normID) {
      setExNorm([String(wave), laserID])
    }
  }, [laserID, normID, setExNorm, wave, updateActiveSpectra])

  const handleNormCheck = (_e, checked) => {
    if (checked) {
      setExNorm([String(wave), laserID])
    } else {
      clearNorm()
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
        <Typography style={{ margin: "8px 10px 3px" }}>Laser</Typography>
        <Box flexGrow={2}>
          <div style={{ margin: "0 12px", minWidth: 200 }}>
            <Typography className={classes.label}>Wavelength</Typography>
            <InputSlider
              value={wave}
              setValue={setWave}
              valueLabelDisplay="auto"
              aria-labelledby="range-slider"
            />
          </div>
        </Box>
        <Box>
          <FormControlLabel
            value="normalize"
            label={<Typography className={classes.label}>Norm em. to this</Typography>}
            labelPlacement="end"
            control={
              <Checkbox
                icon={<RadioButtonUnchecked fontSize="small" />}
                checkedIcon={<RadioButtonChecked fontSize="small" />}
                value="checkedI"
                checked={id.startsWith(normID)}
                onChange={handleNormCheck}
              />
            }
          />
        </Box>
      </Box>
    </div>
  )
})

export default CustomLaserCreator
