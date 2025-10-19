import React from "react"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import Checkbox from "@mui/material/Checkbox"
import FormGroup from "@mui/material/FormGroup"
import FormControlLabel from "@mui/material/FormControlLabel"
import Button from "@mui/material/Button"
import Icon from "@mui/material/Icon"
import SearchIcon from "@mui/icons-material/Search"
import ChartIcon from "@mui/icons-material/InsertChart"
import FileIcon from "@mui/icons-material/GetApp"
import SettingsIcon from "@mui/icons-material/Settings"
import CachedIcon from "@mui/icons-material/Cached"
import ShareIcon from "@mui/icons-material/Share"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"

const useStyles = makeStyles(theme => ({
  root: {
    width: "100%",
    "& .MuiPaper-root ": {
      margin: 48,
      [theme.breakpoints.down("xs")]: {
        margin: 22,
        maxHeight: "calc(100% - 44px)",
      },
    },
    "& h6": {
      color: "#333",
      marginTop: 6,
    },
    "& p": {
      color: "#666",
      marginLeft: "2.2rem",
      marginBottom: 10,
      [theme.breakpoints.down("xs")]: {
        marginLeft: ".4rem",
      },
    },
  },
  headerIcon: {
    marginRight: ".6rem",
    color: "#999",
  },
  button: {
    height: 40,
    marginRight: 20,
  },
  title: {
    color: "black",
    marginBottom: 15,
  },
  footer: {
    marginTop: 15,
  },
  docsButton: {
    textAlign: "center",
    margin: "18px auto 12px",
    "&:hover": {
      color: "inherit",
    },
  },
}))

const WelcomeModal = React.memo(function WelcomeModal({
  open,
  close,
  isNew = false,
  ownerInfo,
}) {
  const classes = useStyles()

  const counter =
    ownerInfo &&
    Object.values(ownerInfo).reduce((acc, next) => {
      acc[next.category] = acc[next.category] + 1 || 1
      return acc
    }, {})
  return (
    <Dialog
      open={open}
      onClose={close}
      aria-labelledby="scroll-dialog-title"
      className={classes.root}
      maxWidth="md"
    >
      <DialogTitle id="scroll-dialog-title" style={{ textAlign: "center" }}>
        Welcome to the
        {isNew ? " new " : ""}
        FPbase Spectra Viewer!
        {Object.values(counter).length > 0 && (
          <div className="stats-list">
            {[
              `${counter.P} proteins`,
              `${counter.D} dyes`,
              `${counter.F} filters`,
              `${counter.L} light sources`,
              `${counter.C} detectors`,
            ].join(" ◦ ")}
          </div>
        )}
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
          &nbsp;to quickly lookup and load any spectrum in the database. Or,
          load any optical config from any
          {" "}
          <a
            href="https://www.fpbase.org/microscopes"
            target="_blank"
            rel="noopener noreferrer"
          >
            FPbase microscope
          </a>
          {" "}
          (including all major
          <strong> filter sets </strong>
          from Chroma, Semrock, Omega, and Zeiss). Try it now!
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
            <ChartIcon />
          </Icon>
          Improved Charts
        </Typography>
        <Typography variant="body1" gutterBottom>
          The chart is faster, less buggy, and handles larger numbers of
          spectra.
        </Typography>
        <Typography variant="h6" gutterBottom>
          <Icon className={classes.headerIcon}>
            <SettingsIcon />
          </Icon>
          Configurable
        </Typography>
        <Typography variant="body1" gutterBottom>
          Change the look and feel with a variety of options in the settings
          menu at the bottom left. Set the X-axis range by clicking and dragging
          or directly enter the max and min values into the inputs. Once zoomed,
          shift-click &amp; drag to pan. See
          {" "}
          <a href="https://help.fpbase.org/tools/spectra-viewer#keyboard-shortcuts">
            documentation
          </a>
          {" "}
          for all keyboard shortcuts.
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
        <div style={{ textAlign: "center" }}>
          <Button
            variant="contained"
            color="inherit"
            href="https://help.fpbase.org/tools/spectra-viewer"
            className={classes.docsButton}
          >
            full documentation
          </Button>
        </div>
      </DialogContent>
      <DialogActions>
        <Button style={{ marginRight: 10 }} onClick={close} color="primary">
          Close
        </Button>
      </DialogActions>
    </Dialog>
  )
})

export default WelcomeModal
