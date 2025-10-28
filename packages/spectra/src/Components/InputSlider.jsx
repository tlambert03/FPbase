import Grid from '@mui/material/Grid'
import Input from '@mui/material/Input'
import Slider from '@mui/material/Slider'
import { makeStyles } from '@mui/styles'

export const useStyles = makeStyles({
  root: {
    width: '100%',
  },
  input: {
    width: 42,
    position: 'relative',
    top: -18,
  },
})

const InputSlider = ({ value, setValue, min = 300, max = 999, step = 1 }) => {
  const classes = useStyles()
  const handleSliderChange = (_event, newValue) => {
    setValue(newValue)
  }
  const handleInputChange = (event) => {
    setValue(event.target.value === '' ? '' : Number(event.target.value))
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
              type: 'number',
              'aria-labelledby': 'input-slider',
            }}
          />
        </Grid>
      </Grid>
    </div>
  )
}

export default InputSlider
