import React, { useState, useReducer, useEffect } from "react";
import ReactDOM from "react-dom";

import SpectraSelectForm from "./SpectraSelect";
import { SpectraViewer, SPECTRUM_CHARTOPTS } from "./ChartOptions";
import client from "./client";
import { FilterContext } from "./context";
import { SPECTRA_LIST } from "./queries";
import { getStorageWithExpire } from "./util";

// state = {
//   spectra: ['list', 'of', 'spectraIDs'],
//   inverted: ['list', 'of', 'spectraIDs'],
//   ecNormed: false,
//   qyNormed: false,
//   logScale: false,
//   extremes: null || (300, 800),
// }

const initialState = { spectrumGroups: {}, series: [] };

function reducer(state, action) {
  switch (action.type) {
    case "updateSpectrumGroups": {
      const newState = { ...state, spectrumGroups: action.payload };
      window.SP = action.payload;
      return newState;
    }
    case "updateSeries": {
      // expects an action.add = [] and action.remove = []
      const newState = { ...state };
      if (action.add && action.add.length) newState.series.push(...action.add);
      if (action.remove && action.remove.length) {
        newState.series = newState.series.filter(item => {
          return !action.remove.includes(item.id);
        });
      }
      window.SERIES = newState.series;
      return newState;
    }
    default:
      throw new Error("Unrecognized reducer action type");
  }
}

const fetchSpectrumGroups = dispatch => {
  const cacheKey = "_spectrumGroups";
  let groups = getStorageWithExpire(cacheKey, 10 * 60);

  if (groups) {
    dispatch({ type: "updateSpectrumGroups", payload: JSON.parse(groups) });
    return;
  }
  client.query({ query: SPECTRA_LIST }).then(({ data }) => {
    if (data) {
      groups = data.spectra.reduce((prev, cur) => {
        const next = prev;
        if (!(cur.owner.slug in next)) {
          next[cur.owner.slug] = {
            label: cur.owner.name,
            category: cur.category,
            spectra: []
          };
        }
        next[cur.owner.slug].spectra.push({
          id: cur.id,
          subtype: cur.subtype
        });
        return next;
      }, {});

      localStorage.setItem(cacheKey, JSON.stringify(groups));
      localStorage.setItem(`${cacheKey}:ts`, Date.now());

      dispatch({ type: "updateSpectrumGroups", payload: groups });
    }
  });
};

const App = () => {
  const [chartSeries, setChartSeries] = useState(SPECTRUM_CHARTOPTS);
  const [state, dispatch] = useReducer(reducer, initialState);

  useEffect(() => {
    fetchSpectrumGroups(dispatch);
  }, []);

  useEffect(() => {
    setChartSeries(state.series);
  }, [state.series]);

  const group = state.spectrumGroups;
  let options = [];
  if (group) {
    options = Object.keys(group)
      .filter(slug => group[slug].category === "D")
      .map(slug => ({
        value: slug,
        label: group[slug].label,
        spectra: group[slug].spectra
      }));
  }

  return (
    <FilterContext.Provider value={{ state, dispatch }}>
      <SpectraViewer series={chartSeries} />
      <SpectraSelectForm options={options} />
    </FilterContext.Provider>
  );
};

const initReactSpectra = elem => {
  ReactDOM.render(<App />, document.getElementById(elem));
};

export default initReactSpectra;
