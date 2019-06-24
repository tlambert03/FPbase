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
import SettingsIcon from "@material-ui/icons/Settings"
import CachedIcon from "@material-ui/icons/Cached"
import ShareIcon from "@material-ui/icons/Share"
import Dialog from "@material-ui/core/Dialog"
import DialogActions from "@material-ui/core/DialogActions"
import DialogContent from "@material-ui/core/DialogContent"
import DialogTitle from "@material-ui/core/DialogTitle"

const useStyles = makeStyles(theme => ({
  root: {
    width: "100%",
    "& .MuiPaper-root ": {
      margin: 48,
      [theme.breakpoints.down("xs")]: {
        margin: 22,
        maxHeight: "calc(100% - 44px)"
      }
    },
    "& h6": {
      color: "#333",
      marginTop: 4
    },
    "& p": {
      color: "#666",
      marginLeft: "2.2rem",
      [theme.breakpoints.down("xs")]: {
        marginLeft: ".4rem"
      }
    }
  },
  headerIcon: {
    marginRight: ".6rem",
    color: "#999"
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

const WelcomeModal = ({ open, checked, handleChange, close, isNew }) => {
  const classes = useStyles()

  return (
    <Dialog
      open={open}
      onClose={close}
      aria-labelledby="scroll-dialog-title"
      className={classes.root}
      maxWidth={"md"}
    >
      <DialogTitle id="scroll-dialog-title">
        {" "}
        Welcome to the {isNew ? "new" : ""} FPbase Spectra Viewer!
      </DialogTitle>
      <DialogContent dividers>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <SearchIcon />
          </Icon>
          Quick Entry
        </Typography>
        <Typography variant="body1" gutterBottom>
          Hit&nbsp;
          <span className="kbd">spacebar</span>
          &nbsp;to quickly lookup and load any spectrum group in the database.
          Try it now!
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <FileIcon />
          </Icon>
          Export &amp; Sharing
        </Typography>
        <Typography variant="body1" gutterBottom>
          The chart may now be exported as PNG, PDF, or SVG vector graphics
          format. Or download all of the current data as CSV format. You may
          also save your chart as a URL for sharing or later recall. All in the
          the share icon (&nbsp;
          <ShareIcon />
          &nbsp;) at the bottom right.
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <SettingsIcon />
          </Icon>
          Configurable
        </Typography>
        <Typography variant="body1" gutterBottom>
          Change the look and feel with a variety of options in the settings
          menu at the bottom left. Set the X range by clicking and dragging,
          shift-click-drag to pan, or directly enter the max and min values.
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <CachedIcon />
          </Icon>
          State Recovery
        </Typography>
        <Typography variant="body1" gutterBottom>
          The state of the viewer in any tab will persist across browser refresh
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
            style={{ paddingTop: 8 }}
            control={
              <Checkbox
                checked={checked}
                color="primary"
                onChange={handleChange}
                value="checked"
              />
            }
            label="Don't show on load"
          />
        </FormGroup>
        <Button onClick={close} color="primary">
          Close
        </Button>
      </DialogActions>
    </Dialog>
  )
}

WelcomeModal.defaultProps = {
  isNew: false
}
export default WelcomeModal
