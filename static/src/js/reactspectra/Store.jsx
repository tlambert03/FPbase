import React, { useReducer, useEffect } from "react";
import qs from "qs";
import client from "./client";
import { SPECTRA_LIST } from "./queries";
import { getStorageWithExpire, DEFAULT_EXPIRY } from "./util";

const initialState = {
  owners: {},
  ownerCategories: {},
  spectraInfo: {},
  currentSpectra: [],
  inverted: [],
  ecNormed: false,
  qyNormed: false,
  loading: true,
  tab: 0
};

const AppContext = React.createContext(initialState);

function reducer(state, action) {
  switch (action.type) {
    case "UPDATE_SPECTRA": {
      let { currentSpectra } = state;
      if (action.remove) {
        currentSpectra = currentSpectra.filter(i => !action.remove.includes(i));
      }
      if (action.add) {
        currentSpectra = [...currentSpectra, ...action.add];
      }
      return { ...state, currentSpectra };
    }
    case "REMOVE_SPECTRA": {
      const spectra = state.currentSpectra.filter(
        i => !action.payload.includes(i)
      );
      return { ...state, currentSpectra: spectra };
    }
    case "ADD_SPECTRA": {
      return {
        ...state,
        currentSpectra: [...state.currentSpectra, ...action.payload]
      };
    }
    case "UPDATE": {
      const newState = { ...state };
      Object.keys(action).forEach(prop => {
        if (prop !== "type" && prop in newState) {
          newState[prop] = action[prop];
        }
      });
      window.STATE = newState;
      return newState;
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

const initialize = async dispatch => {
  // look for params in the URL
  const urlParams = qs.parse(window.location.search.substr(1));
  const urlSeries = urlParams.series || urlParams.s;
  if (urlSeries) {
    const ids = [...new Set(urlSeries.split(","))];
    dispatch({ type: "ADD_SPECTRA", payload: ids });
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
            spectra: []
          };
        }
        prev.owners[cur.owner.slug].spectra.push({
          id: cur.id,
          subtype: cur.subtype
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
  window.owners = owners;
  window.spectraInfo = spectraInfo;
  dispatch({
    type: "UPDATE",
    owners,
    spectraInfo,
    ownerCategories: seperateOwnerCategories(owners),
    loading: false
  });
};

const Store = ({ children }) => {
  const [state, dispatch] = useReducer(reducer, initialState);

  useEffect(() => {
    initialize(dispatch);
  }, []);

  return (
    <AppContext.Provider value={{ state, dispatch }}>
      {children}
    </AppContext.Provider>
  );
};

export { Store, AppContext };
