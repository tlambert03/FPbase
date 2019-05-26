import React, { useState, useContext, useEffect } from "react";
import Select from "react-select";
import Box from "@material-ui/core/Box";
import ToggleButton from "@material-ui/lab/ToggleButton";
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup";
import { makeStyles } from "@material-ui/core/styles";
import { AppContext } from "./Store";
import { ID } from "./util";

const useStyles = makeStyles(theme => ({
  toggleButton: { height: "38px" },
  toggleButtonGroup: { marginLeft: "5px" }
}));

export const customFilterOption = (option, rawInput) => {
  const words = rawInput.split(" ");
  return words.reduce(
    (acc, cur) => acc && option.label.toLowerCase().includes(cur.toLowerCase()),
    true
  );
};

function subtypeSorter(a, b) {
  const TYPE_ORDER = ["AB", "A_2P", "EX", "EM"];
  const upperA = a.subtype.toUpperCase();
  const upperB = b.subtype.toUpperCase();
  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1;
  return -1;
}

// where value = {value: "ownerslug", lable: "ownername", spectra: [Array, of, IDs]}
const SpectrumSelector = ({
  options,
  initialValue,
  initialSubtypes,
  initialAvailableSubtypes
}) => {
  const [value, setValue] = useState(initialValue);
  const [subtypes, setSubtypes] = useState(initialSubtypes);
  const [availableSubtypes, setAvailableSubtypes] = useState(
    initialAvailableSubtypes
  );
  const { dispatch } = useContext(AppContext);

  // when the spectrum selector changes
  const handleOwnerChange = e => {
    // if it's the same as the previous value do nothing
    if (e && e.value === value) return;
    dispatch({
      type: "UPDATE_SPECTRA",
      add: e && e.spectra.map(({ id }) => id),
      remove: value && value.spectra.map(({ id }) => id)
    });
    setSubtypes((e && e.spectra) || []);
    setValue(e);
  };

  // clean up on unmount
  useEffect(() => {
    setAvailableSubtypes((value && value.spectra) || []);
    return () => {
      if (value)
        dispatch({
          type: "REMOVE_SPECTRA",
          payload: value.spectra.map(({ id }) => id)
        });
    };
  }, [dispatch, value]);

  const handleSelectionChange = (e, newSelection) => {
    setSubtypes(
      availableSubtypes.filter(({ id }) => newSelection.includes(id))
    );
  };

  return (
    <Box display="flex">
      <Box flexGrow={1}>
        <Select
          isClearable
          defaultValue={value}
          placeholder="Type to search..."
          filterOption={customFilterOption}
          onChange={handleOwnerChange}
          options={options}
        />
      </Box>
      {availableSubtypes.length > 1 && (
        <SubtypeSelector
          subtypes={availableSubtypes}
          onChange={handleSelectionChange}
          intialValue={subtypes.map(({ id }) => id)}
        />
      )}
    </Box>
  );
};

SpectrumSelector.defaultProps = {
  defaultSubtypes: []
};

const SubtypeSelector = ({ subtypes, onChange, intialValue }) => {
  const { dispatch } = useContext(AppContext);
  const classes = useStyles();

  const handleClick = e => {
    const elem = e.target.closest("button");
    const checked = !elem.classList.contains("Mui-selected");
    const type = checked ? "ADD_SPECTRA" : "REMOVE_SPECTRA";
    dispatch({ type, payload: [elem.value] });
  };

  subtypes.sort(subtypeSorter);
  return (
    <Box>
      <ToggleButtonGroup
        value={intialValue}
        onChange={onChange}
        size="small"
        className={classes.toggleButtonGroup}
      >
        {subtypes.map(st => (
          <ToggleButton
            key={st.id}
            value={st.id}
            onClick={handleClick}
            className={classes.toggleButton}
          >
            {st.subtype.replace(/^A_/g, "")}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
    </Box>
  );
};

export default SpectrumSelector;
