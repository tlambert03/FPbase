import React, { useState, useEffect } from "react";
import Modal from "@material-ui/core/Modal";
import Typography from "@material-ui/core/Typography";
import { makeStyles } from "@material-ui/core/styles";
import Checkbox from "@material-ui/core/Checkbox";
import FormGroup from "@material-ui/core/FormGroup";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import Button from "@material-ui/core/Button";
import Box from "@material-ui/core/Box";

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
  paper: {
    position: "absolute",
    width: "80%",
    backgroundColor: theme.palette.background.paper,
    boxShadow: theme.shadows[5],
    padding: 30,
    paddingTop: 25,
    outline: "none"
  },
  button: {
    height: 40,
    marginRight: 20
  },
  title: {
    color: "black",
    marginBottom: 15
  }
}));

const SearchModal = () => {
  const [modalStyle] = useState(getModalStyle);
  const [checked, setChecked] = useState(false);
  const classes = useStyles();
  const [open, setOpen] = useState(false);

  const storageKey = "_hideFPbaseSpectraWelcome";
  useEffect(() => {
    const show = localStorage.getItem(storageKey);
    if (show !== "true") {
      setOpen(true);
    }
  }, []);

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
      <div style={modalStyle} className={classes.paper}>
        <Typography
          variant="h6"
          id="welcome-modal-title"
          className={classes.title}
        >
          Welcome to the new FPbase Spectra Viewer!
        </Typography>
        <Box display="flex">
          <Button
            variant="contained"
            color="primary"
            className={classes.button}
            onClick={() => setOpen(false)}
          >
            Got it
          </Button>
          <FormGroup row>
            <FormControlLabel
              // prettier-ignore
              control={<Checkbox checked={checked} color="primary" onChange={handleChange} value="checked" />}
              label="Don't show this again"
            />
          </FormGroup>
        </Box>
      </div>
    </Modal>
  );
};

export default SearchModal;
