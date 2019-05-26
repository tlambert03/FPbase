import React, { useState, useContext } from "react";
import Select from "react-select";
import Box from "@material-ui/core/Box";
import IconButton from "@material-ui/core/IconButton";
import Button from "@material-ui/core/Button";
import DeleteIcon from "@material-ui/icons/Delete";
import AddIcon from "@material-ui/icons/Add";
import { makeStyles } from "@material-ui/core/styles";
import { ID } from "./util";
import SpectrumSelector from "./SpectrumSelector";

const useStyles = makeStyles(theme => ({
  button: {
    "&:hover": {
      // you want this to be the same as the backgroundColor above
      backgroundColor: "#f8f8f8"
    }
  }
}));

const SpectrumSelectGroup = ({ options, hint, formState }) => {
  const [selectors, setSelectors] = useState(formState);
  const classes = useStyles();

  const addRow = () => setSelectors([...selectors, { id: ID() }]);
  const removeSelector = selID => {
    setSelectors(selectors.filter(({ id }) => id !== selID));
  };

  const makeSelector = selector => {
    let initial = null;
    if (selector.value) {
      initial = options.find(({ value }) => value === selector.value);
    }
    return (
      <SpectrumSelector
        id={selector.id}
        options={options}
        initialValue={initial}
        initialSubtypes={selector.spectra}
        initialAvailableSubtypes={(initial && initial.spectra) || []}
      />
    );
  };

  return (
    <>
      {options.length > 0 ? (
        selectors.map(selector => (
          <div style={{ width: "100%" }} key={selector.id}>
            <Box display="flex">
              <Box flexGrow={1} style={{ paddingTop: "6px" }}>
                {makeSelector(selector)}
              </Box>
              {selectors.length > 1 ? (
                <Box>
                  <IconButton
                    aria-label="Delete"
                    color="secondary"
                    className={classes.button}
                    onClick={() => removeSelector(selector.id)}
                  >
                    <DeleteIcon />
                  </IconButton>
                </Box>
              ) : (
                <></>
              )}
            </Box>
          </div>
        ))
      ) : (
        <Select placeholder="Loading..." isLoading />
      )}
      <Button
        variant="contained"
        color="primary"
        className="mt-2"
        onClick={addRow}
      >
        {`Add ${hint || "item"}`}
        <AddIcon />
      </Button>
    </>
  );
};

SpectrumSelectGroup.defaultProps = {
  // one id: is required otherwise there will be no rows
  formState: [{ id: ID(), value: null, subtypes: [] }]
};

export default SpectrumSelectGroup;
