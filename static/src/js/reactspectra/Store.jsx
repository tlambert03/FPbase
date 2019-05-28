import React, { useReducer, useEffect } from "react";
import PropTypes from "prop-types";
import qs from "qs";
import client from "./client";
import { SPECTRA_LIST } from "./queries";
import {
  getStorageWithExpire,
  DEFAULT_EXPIRY,
  ID,
  emptyFormSelector
} from "./util";

const initialState = {
  owners: {}, // mostly duplicate of ownerCategories, ownerSlug lookup
  ownerCategories: {}, // all possible owners, by organized by spectrum Category
  spectraInfo: {},
  formState: {},
  currentSpectra: [],
  inverted: [],
  ecNormed: false,
  qyNormed: false,
  loading: true,
  pushState: true,
  tab: 0
};

const AppContext = React.createContext(initialState);

// This function is responsible for making sure the form (minimally)
// reflects the spectra ids in State.currentSpectra
// it does not prevent the form from having extra rows... or duplicate owners
const stateToSpectraForms = (
  currentSpectra,
  spectraInfo,
  formState,
  owners
) => {
  // get spectra info for each item in the currentSpectra
  const info = currentSpectra
    // .filter(id => {
    //   if (Object.prototype.hasOwnProperty.call(spectraInfo, id)) {
    //     return true;
    //   }
    //   console.error(`could not get info for spectrum ID ${id}`) // eslint-disable-line
    //   return false;
    // })
    .map(id => spectraInfo[id]);

  // figure out what owner and category they all belong to
  const activeOwners = info.reduce((acc, obj) => {
    acc[obj.category] = (acc[obj.category] || new Set()).add(obj.owner);
    return acc;
  }, {});

  // building a new formState
  const newForm = { ...formState };

  // make sure all of the active owners are somewhere on the form
  Object.keys(activeOwners).forEach(key => {
    activeOwners[key].forEach(owner => {
      newForm[key] = newForm[key] || [];
      if (!newForm[key].some(({ value }) => value === owner)) {
        newForm[key] = [...newForm[key], { id: ID(), ...owners[owner] }];
      }
    });
  });

  // iterate through the existing form
  Object.keys(newForm).forEach(key => {
    // then make sure all of the currentSpectra are marked as active on the form
    newForm[key].forEach((elem, idx) => {
      if (elem.spectra) {
        newForm[key][idx].spectra = elem.spectra.map(i => {
          // eslint-disable-next-line no-param-reassign
          i.active = currentSpectra.includes(i.id);
          return i;
        });
      }
    });
    // put empty search fields at the end
    newForm[key].sort((a, b) => (a.value && !b.value ? -1 : 0));
  });
  ["P", "D", "L", "F", "C"].forEach(letter => {
    if (!(letter in newForm)) {
      newForm[letter] = [emptyFormSelector()];
    }
  });
  return newForm;
};

const addAndRemoveFromArray = (array, action) => {
  let newArray = array;
  if (action.remove) {
    newArray = newArray.filter(i => !action.remove.includes(i));
  }
  if (action.add) {
    newArray = [...newArray, ...action.add];
  }
  return [...new Set(newArray)];
};

function reducer(state, action) {
  switch (action.type) {
    case "UPDATE_SPECTRA": {
      const currentSpectra = addAndRemoveFromArray(
        state.currentSpectra,
        action
      );
      return { ...state, currentSpectra };
    }
    case "UPDATE": {
      const newState = { ...state };
      Object.keys(action).forEach(prop => {
        if (prop !== "type" && prop in newState) {
          newState[prop] = action[prop];
        }
      });
      return newState;
    }
    case "CHANGE_TAB": {
      const newState = { ...state };
      newState.tab = action.payload;
      return newState;
    }
    case "ADD_FORM_ROW": {
      const newState = { ...state };
      const { category } = action;
      const currentRow = newState.formState[category] || [];
      newState.formState[category] = [...currentRow, emptyFormSelector()];
      return newState;
    }
    case "REMOVE_FORM_ROW": {
      const newState = { ...state };
      const { category, id } = action;
      const currentRow = newState.formState[category];
      newState.formState[category] = currentRow.filter(item => item.id !== id);
      return newState;
    }
    case "CHANGE_FORM_OWNER": {
      // we change both the spectra and the owner here
      // so that we don't get a duplicated form element
      const newState = { ...state };
      const { category, id, newValue } = action;
      const currentRow = newState.formState[category];
      const idx = currentRow.findIndex(i => i.id === id);
      const newOwner = state.owners[newValue];
      const oldOwner = newState.formState[category][idx];
      const changeAction = {
        add: newOwner && newOwner.spectra && newOwner.spectra.map(i => i.id),
        remove: oldOwner && oldOwner.spectra && oldOwner.spectra.map(i => i.id)
      };
      // change the form
      if (idx > -1) {
        newState.formState[category][idx] = {
          id,
          ...newOwner
        };
      }
      // change the currentSpectra
      const currentSpectra = addAndRemoveFromArray(
        state.currentSpectra,
        changeAction
      );
      return { ...newState, currentSpectra };
    }
    case "UPDATE_FORM": {
      const { currentSpectra, spectraInfo, formState, owners } = state;
      const newState = { ...state };
      newState.formState = stateToSpectraForms(
        currentSpectra,
        spectraInfo,
        formState,
        owners
      );
      return newState;
    }
    case "RESET": {
      return { ...initialState, formState: {} };
    }
    default:
      throw new Error(`Unrecognized reducer action type: ${action.type}`);
  }
}

const seperateOwnerCategories = owners => {
  return Object.keys(owners).reduce((accum, slug) => {
    const { category, label, spectra } = owners[slug];
    if (!Object.prototype.hasOwnProperty.call(accum, category)) {
      // eslint-disable-next-line no-param-reassign
      accum[category] = [];
    }
    accum[category].push({
      value: slug,
      label,
      spectra
    });
    return accum;
  }, {});
};

const initialize = async (dispatch, parseURL = true) => {
  // look for params in the URL
  if (parseURL) {
    const urlParams = qs.parse(window.location.search.substr(1));
    const urlSeries = urlParams.series || urlParams.s;
    if (urlSeries) {
      const ids = [...new Set(urlSeries.split(","))];
      dispatch({ type: "UPDATE_SPECTRA", add: ids });
    }
    const tab = urlParams.tab || urlParams.t;
    if (tab) {
      dispatch({ type: "UPDATE", tab: +tab });
    }
  }

  // Grab available spectra ids and owner slugs from storage or server
  const cacheKey = "_availableSpectra";
  let cachedData = getStorageWithExpire(cacheKey, DEFAULT_EXPIRY);
  if (cachedData) {
    cachedData = JSON.parse(cachedData);
  } else {
    const { data } = await client.query({ query: SPECTRA_LIST });
    cachedData = data.spectra.reduce(
      (prev, cur) => {
        if (
          !Object.prototype.hasOwnProperty.call(prev.owners, cur.owner.slug)
        ) {
          // eslint-disable-next-line no-param-reassign
          prev.owners[cur.owner.slug] = {
            category: cur.category,
            label: cur.owner.name,
            spectra: [],
            value: cur.owner.slug,
            url: cur.owner.url
          };
        }
        prev.owners[cur.owner.slug].spectra.push({
          id: cur.id,
          subtype: cur.subtype,
          active: true
        });
        // eslint-disable-next-line no-param-reassign
        prev.spectraInfo[cur.id] = {
          subtype: cur.subtype,
          owner: cur.owner.slug,
          label: cur.owner.name,
          category: cur.category
        };
        return prev;
      },
      { owners: {}, spectraInfo: {} }
    );
    localStorage.setItem(cacheKey, JSON.stringify(cachedData));
    localStorage.setItem(`${cacheKey}:ts`, Date.now());
  }
  const { owners, spectraInfo } = cachedData;

  dispatch({
    type: "UPDATE",
    owners,
    spectraInfo,
    ownerCategories: seperateOwnerCategories(owners)
  });
  // next line prevents an empty form item at the beginning
  // before all of the populated form items
  dispatch({ type: "UPDATE_FORM" });
  dispatch({ type: "UPDATE", loading: false });
};

const Store = ({ children }) => {
  const [state, dispatch] = useReducer(reducer, initialState);

  useEffect(() => {
    initialize(dispatch);
  }, []);

  // This effect keeps the form consistent with the
  // spectra in state.currentSpectra
  useEffect(() => {
    if (!state.loading) dispatch({ type: "UPDATE_FORM" });
  }, [state.currentSpectra]) // eslint-disable-line

  return (
    <AppContext.Provider value={{ state, dispatch }}>
      {children}
    </AppContext.Provider>
  );
};

Store.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
};

export { Store, AppContext, initialize };
