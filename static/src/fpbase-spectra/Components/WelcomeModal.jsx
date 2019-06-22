import React, { useState } from "react"
import Typography from "@material-ui/core/Typography"
import { makeStyles } from "@material-ui/core/styles"
import Checkbox from "@material-ui/core/Checkbox"
import FormGroup from "@material-ui/core/FormGroup"
import FormControlLabel from "@material-ui/core/FormControlLabel"
import Button from "@material-ui/core/Button"
import Icon from "@material-ui/core/Icon"
import SearchIcon from "@material-ui/icons/Search"
import ChartIcon from "@material-ui/icons/InsertChart"
import FileIcon from "@material-ui/icons/GetApp"
import CachedIcon from "@material-ui/icons/Cached"

import Dialog from "@material-ui/core/Dialog"
import DialogActions from "@material-ui/core/DialogActions"
import DialogContent from "@material-ui/core/DialogContent"
import DialogTitle from "@material-ui/core/DialogTitle"

const useStyles = makeStyles(theme => ({
  root: {

    "& h6": {
      color: "#333"
    },
    "& p": {
      color: "#666",
      marginLeft: "2.2rem"
    }
  },
  headerIcon: {
    marginRight: ".6rem",
    color: "#999",
    top: 5,
    position: "relative"
  },
  button: {
    height: 40,
    marginRight: 20
  },
  title: {
    color: "black",
    marginBottom: 15
  },
  footer: {
    marginTop: 15
  }
}))

const SearchModal = () => {
  const classes = useStyles()
  const storageKey = "_hideFPbaseSpectraWelcome"
  const hide = localStorage.getItem(storageKey) === "true"
  const [checked, setChecked] = useState(hide)
  const [open, setOpen] = useState(!hide)

  const handleChange = e => {
    localStorage.setItem(storageKey, e.target.checked)
    setChecked(e.target.checked)
  }

  return (
    <Dialog
      open={open}
      onClose={() => setOpen(false)}
      aria-labelledby="scroll-dialog-title"
      className={classes.root}
      fullWidth
      maxWidth={'md'}
    >
      <DialogTitle id="scroll-dialog-title">
        {" "}
        Welcome to the new FPbase Spectra Viewer!
      </DialogTitle>
      <DialogContent dividers>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <SearchIcon />
          </Icon>
          Quick Entry
        </Typography>
        <Typography variant="body1" gutterBottom>
          Hit the
          <span className="kbd">L</span>
          &nbsp;key to quickly lookup and load any spectrum group in the
          database. Try it now!
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <FileIcon />
          </Icon>
          Export Options
        </Typography>
        <Typography variant="body1" gutterBottom>
          The chart may now be exported as PNG, PDF, or SVG vector graphics
          format. You can also download all of the data in the current chart in
          CSV format. (Use the
          <Icon
            className="fas fa-bars"
            style={{ fontSize: "0.93rem", margin: "0 5px" }}
          />
          context menu at the top right)
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <CachedIcon />
          </Icon>
          State Recovery
        </Typography>
        <Typography variant="body1" gutterBottom>
          The state of the spectra viewer is captured in the URL. So you can
          bookmark or refresh without losing your work. Or share a specific
          arrangement of spectra with someone else
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <ChartIcon />
          </Icon>
          Better Charts
        </Typography>
        <Typography variant="body1" gutterBottom>
          The chart is faster, less buggy, and handles larger numbers of
          spectra.
        </Typography>
      </DialogContent>
      <DialogActions>
        <FormGroup row>
          <FormControlLabel
            control={
              <Checkbox
                checked={checked}
                color="primary"
                onChange={handleChange}
                value="checked"
              />
            }
            label="Don't show again"
          />
        </FormGroup>
        <Button onClick={() => setOpen(false)} color="primary">
          Close
        </Button>
      </DialogActions>
    </Dialog>
  )
}

export default SearchModal
