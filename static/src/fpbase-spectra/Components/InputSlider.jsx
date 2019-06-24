import React from "react"
import Slider from "@material-ui/lab/Slider"
import Grid from "@material-ui/core/Grid"
import Input from "@material-ui/core/Input"
import { useStyles } from "./CustomFilterCreator"

const InputSlider = ({ value, setValue, min, max, step }) => {
  const classes = useStyles()
  const handleSliderChange = (event, newValue) => {
    setValue(newValue)
  }
  const handleInputChange = event => {
    setValue(event.target.value === "" ? "" : Number(event.target.value))
  }
  const handleBlur = () => {
    if (value < min) {
      setValue(min)
    } else if (value > max) {
      setValue(max)
    }
  }
  return (
    <div className={classes.root}>
      <Grid container spacing={2} alignItems="center">
        <Grid item xs>
          <Slider
            value={+value}
            onChange={handleSliderChange}
            aria-labelledby="input-slider"
            min={min}
            max={max}
            step={step}
          />
        </Grid>
        <Grid item>
          <Input
            className={classes.input}
            value={value}
            margin="dense"
            onChange={handleInputChange}
            onBlur={handleBlur}
            inputProps={{
              step,
              min,
              max,
              type: "number",
              "aria-labelledby": "input-slider"
            }}
          />
        </Grid>
      </Grid>
    </div>
  )
}

InputSlider.defaultProps = {
  min: 300,
  max: 1000,
  step: 1
}

export default InputSlider
