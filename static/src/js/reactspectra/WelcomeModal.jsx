import React, { useState } from "react";
import Modal from "@material-ui/core/Modal";
import Typography from "@material-ui/core/Typography";
import { makeStyles } from "@material-ui/core/styles";
import Checkbox from "@material-ui/core/Checkbox";
import FormGroup from "@material-ui/core/FormGroup";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import Button from "@material-ui/core/Button";
import Box from "@material-ui/core/Box";
import Divider from "@material-ui/core/Divider";
import Icon from "@material-ui/core/Icon";
import SearchIcon from "@material-ui/icons/Search";
import ChartIcon from "@material-ui/icons/InsertChart";
import FileIcon from "@material-ui/icons/GetApp";
import CachedIcon from "@material-ui/icons/Cached";

function getModalStyle() {
  const top = 40;
  const left = 42;

  return {
    top: `${top}%`,
    left: `${left}%`,
    transform: `translate(-${top}%, -${left}%)`
  };
}

const useStyles = makeStyles(theme => ({
  root: {
    width: "100%",
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
    color: "#999"
  },
  paper: {
    position: "absolute",
    width: "65%",
    backgroundColor: theme.palette.background.paper,
    boxShadow: theme.shadows[5],
    padding: 35,
    paddingTop: 30,
    paddingBottom: 15,
    outline: "none"
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
}));

const SearchModal = () => {
  const [modalStyle] = useState(getModalStyle);
  const classes = useStyles();
  const storageKey = "_hideFPbaseSpectraWelcome";
  const hide = localStorage.getItem(storageKey) === "true";
  const [checked, setChecked] = useState(hide);
  const [open, setOpen] = useState(!hide);

  const handleChange = e => {
    localStorage.setItem(storageKey, e.target.checked);
    setChecked(e.target.checked);
  };

  return (
    <Modal
      aria-labelledby="welcome-modal"
      aria-describedby="welcome-modal"
      open={open}
      onClose={() => setOpen(false)}
    >
      <div className={classes.root}>
        <div style={modalStyle} className={classes.paper}>
          <Typography variant="h5" id="welcome-modal" className={classes.title}>
            Welcome to the new FPbase Spectra Viewer!
          </Typography>

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
            More Exporting Options
          </Typography>
          <Typography variant="body1" gutterBottom>
            The chart may now be exported as PNG, PDF, or SVG vector graphics
            format. You can also download all of the data in the current chart
            in CSV format. (Use the
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
          <Divider style={{ marginTop: 20 }} />
          <Box display="flex" className={classes.footer}>
            <Button
              variant="contained"
              color="default"
              className={classes.button}
              onClick={() => setOpen(false)}
            >
              Close
            </Button>
            <FormGroup row>
              <FormControlLabel
                // prettier-ignore
                control={<Checkbox checked={checked} color="default" onChange={handleChange} value="checked" />}
                label="Don't show this again"
              />
            </FormGroup>
          </Box>
        </div>
      </div>
    </Modal>
  );
};

export default SearchModal;
