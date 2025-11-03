import RadioButtonChecked from "@mui/icons-material/RadioButtonChecked"
import RadioButtonUnchecked from "@mui/icons-material/RadioButtonUnchecked"
import Box from "@mui/material/Box"
import Checkbox from "@mui/material/Checkbox"
import FormControlLabel from "@mui/material/FormControlLabel"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import { memo, useEffect } from "react"
import { useSpectraStore } from "../store/spectraStore"
import InputSlider from "./InputSlider"

const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444",
  },
})

// id is stable like "$cl1"
const CustomLaserCreator = memo(function CustomLaserCreator({ id, normID, setExNorm, clearNorm }) {
  const classes = useStyles()

  // Read params directly from store (single source of truth)
  const params = useSpectraStore((state) => state.customLasers[id])
  const updateCustomLaser = useSpectraStore((state) => state.updateCustomLaser)

  // Default value if params don't exist yet
  const wavelength = params?.wavelength || 488

  // Update exNorm when wavelength changes (if this laser is the norm target)
  useEffect(() => {
    if (id === normID) {
      setExNorm([wavelength, id])
    }
  }, [id, normID, setExNorm, wavelength])

  const handleNormCheck = (_e, checked) => {
    if (checked) {
      setExNorm([wavelength, id])
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
              value={wavelength}
              setValue={(value) => updateCustomLaser(id, { wavelength: value })}
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
