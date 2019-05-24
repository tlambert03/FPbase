import React, { useState, useContext, useEffect } from "react";
import Select from "react-select";
import Box from "@material-ui/core/Box";
import ToggleButton from "@material-ui/lab/ToggleButton";
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup";
import IconButton from "@material-ui/core/IconButton";
import Button from "@material-ui/core/Button";
import DeleteIcon from "@material-ui/icons/Delete";
import AddIcon from "@material-ui/icons/Add";
import { makeStyles } from "@material-ui/core/styles";
import { ID, updateSpectra } from "./util";
import { FilterContext } from "./context";

const useStyles = makeStyles(theme => ({
  button: {},
  toggleButton: { height: "38px" },
  toggleButtonGroup: { marginLeft: "5px" },
  selector: {
    paddingTop: "6px"
  }
}));

const customFilterOption = (option, rawInput) => {
  const words = rawInput.split(" ");
  return words.reduce(
    (acc, cur) => acc && option.label.toLowerCase().includes(cur.toLowerCase()),
    true
  );
};

const SpectraSelectForm = ({ options }) => {
  const [selectors, setSelectors] = useState([{ id: ID() }]);
  const classes = useStyles();

  const addRow = () => setSelectors([...selectors, { id: ID() }]);
  const removeSelector = selID => {
    setSelectors(selectors.filter(({ id }) => id !== selID));
  };
  return (
    <>
      {options.length > 0 ? (
        selectors.map(selector => (
          <div style={{ width: "100%" }} key={selector.id}>
            <Box display="flex">
              <Box flexGrow={1} className={classes.selector}>
                <SpectrumSelect id={selector.id} options={options} />
              </Box>
              <Box>
                <IconButton
                  aria-label="Delete"
                  className={classes.button}
                  onClick={() => removeSelector(selector.id)}
                >
                  <DeleteIcon />
                </IconButton>
              </Box>
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
        Add item
        <AddIcon />
      </Button>
    </>
  );
};

function dyeCompare(a, b) {
  const TYPE_ORDER = ["AB", "A_2P", "EX", "EM"];
  const upperA = a.toUpperCase();
  const upperB = b.toUpperCase();

  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1;
  return -1;
}

const SpectrumSelect = ({ options }) => {
  const [value, setValue] = useState();
  const { dispatch } = useContext(FilterContext);
  const [selectedTypes, setSelectedTypes] = React.useState(["ABS", "EX", "EM"]);
  const classes = useStyles();

  const handleType = (event, newTypes) => {
    setSelectedTypes(newTypes);
  };

  // when the selector changes
  const handleChange = e => {
    if (e && e.value === value) return;
    setValue(e || null);
    updateSpectra({
      toAdd: e && e.spectra,
      toRemove: value && value.spectra,
      dispatch
    });
  };

  // clean up on unmount
  useEffect(() => {
    return () => {
      updateSpectra({
        toAdd: null,
        toRemove: value && value.spectra,
        dispatch
      });
    };
  }, [dispatch, value]);

  let subtypes = [];
  if (value && value.spectra && value.spectra.length) {
    subtypes = value.spectra.map(e => e.subtype);
  }
  subtypes.sort(dyeCompare);

  return (
    <>
      <Box display="flex">
        <Box flexGrow={1}>
          <Select
            isClearable
            defaultValue={value}
            placeholder="Type to search..."
            filterOption={customFilterOption}
            onChange={handleChange}
            options={options.map(_opt => {
              const opt = { ..._opt };
              if (opt.value === value) {
                opt.isDisabled = false;
              }
              return opt;
            })}
          />
        </Box>
        {subtypes.length > 1 && (
          <SubtypeSelector
            subtypes={subtypes}
            onChange={handleType}
            value={selectedTypes}
          />
        )}
      </Box>
    </>
  );
};

const SubtypeSelector = ({ subtypes, onChange, value }) => {
  const classes = useStyles();

  return (
    <Box>
      <ToggleButtonGroup
        value={value}
        size="small"
        onChange={onChange}
        className={classes.toggleButtonGroup}
      >
        {subtypes.map(st => (
          <ToggleButton
            key={st}
            value={st}
            color="primary"
            className={classes.toggleButton}
          >
            {st}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
    </Box>
  );
};

export default SpectraSelectForm;
