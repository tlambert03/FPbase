import React, { useReducer, useEffect } from "react";
import qs from "qs";
import client from "./client";
import { SPECTRA_LIST } from "./queries";
import { getStorageWithExpire, fetchSpectraList, DEFAULT_EXPIRY } from "./util";
import { reducer } from "./reducer";

const initialState = {
  owners: {},
  spectra: {},
  series: [],
  inverted: [],
  ecNormed: false,
  qyNormed: false,
  loading: true,
  tab: 0
};

const AppContext = React.createContext(initialState);

const initialize = async dispatch => {
  // look for params in the URL
  const urlParams = qs.parse(window.location.search.substr(1));
  const urlSeries = urlParams.series || urlParams.s;
  if (urlSeries) {
    const ids = [...new Set(urlSeries.split(","))];
    const newSpectra = await fetchSpectraList(ids);
    dispatch({ type: "ADD_SERIES", toAdd: newSpectra });
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
        const next = { ...prev };
        if (!(cur.owner.slug in next.owners)) {
          next.owners[cur.owner.slug] = {
            category: cur.category,
            label: cur.owner.name,
            spectra: []
          };
        }
        next.owners[cur.owner.slug].spectra.push(cur.id);
        next.spectra[cur.id] = {
          subtype: cur.subtype,
          owner: cur.owner.slug,
          category: cur.category
        };
        return next;
      },
      { owners: {}, spectra: {} }
    );
    localStorage.setItem(cacheKey, JSON.stringify(cachedData));
    localStorage.setItem(`${cacheKey}:ts`, Date.now());
  }
  const { owners, spectra } = cachedData;
  dispatch({
    type: "UPDATE",
    owners,
    spectra,
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
