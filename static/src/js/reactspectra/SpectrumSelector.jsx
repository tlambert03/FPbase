import React, { useContext, useEffect } from "react";
import PropTypes from "prop-types";
import Select from "react-select";
import Box from "@material-ui/core/Box";
import ToggleButton from "@material-ui/lab/ToggleButton";
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup";
import { makeStyles } from "@material-ui/core/styles";
import { AppContext } from "./Store";

const useStyles = makeStyles(() => ({
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
const SpectrumSelector = ({ options, selector, category }) => {
  const {
    state: { formState },
    dispatch
  } = useContext(AppContext);

  const value = options.find(opt => opt.value === selector.value);
  const subtypes = (value && value.spectra) || [];

  // when the spectrum selector changes
  const handleOwnerChange = e => {
    // if it's the same as the previous value do nothing
    const newValue = e && e.value;
    if (newValue === value) return;
    dispatch({
      type: "CHANGE_FORM_OWNER",
      id: selector.id,
      category,
      newValue
    });
  };

  // clean up on unmount
  useEffect(() => {
    return () => {
      if (value)
        dispatch({
          type: "UPDATE_SPECTRA",
          remove: value.spectra.map(({ id }) => id)
        });
    };
  }, [dispatch, value]);

  const otherOwners = formState[category]
    .filter(i => i.value && i.value !== selector.value)
    .map(i => i.value);
  const myOptions = options.filter(opt => !otherOwners.includes(opt.value));

  return (
    <Box display="flex">
      <Box flexGrow={1}>
        <Select
          isClearable
          defaultValue={value}
          placeholder="Type to search..."
          filterOption={customFilterOption}
          onChange={handleOwnerChange}
          options={myOptions}
        />
      </Box>
      {subtypes.length > 1 && <SubtypeSelector subtypes={subtypes} />}
    </Box>
  );
};

SpectrumSelector.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired,
  selector: PropTypes.arrayOf(PropTypes.object).isRequired,
  category: PropTypes.string.isRequired
};

const SubtypeSelector = ({ subtypes }) => {
  const { dispatch } = useContext(AppContext);
  const classes = useStyles();

  const handleClick = e => {
    const elem = e.target.closest("button");
    const checked = !elem.classList.contains("Mui-selected");
    const action = { type: "UPDATE_SPECTRA" };
    action[checked ? "add" : "remove"] = [elem.value];
    dispatch(action);
  };

  subtypes.sort(subtypeSorter);
  return (
    <Box>
      <ToggleButtonGroup
        value={subtypes.filter(i => i.active).map(i => i.id)}
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

SubtypeSelector.propTypes = {
  subtypes: PropTypes.arrayOf(PropTypes.objectOf(PropTypes.string)).isRequired
};

export default SpectrumSelector;
