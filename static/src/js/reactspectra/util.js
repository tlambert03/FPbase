import { GET_SPECTRUM } from "./queries";
import client from "./client";

const ID = () => {
  // Math.random should be unique because of its seeding algorithm.
  // Convert it to base 36 (numbers + letters), and grab the first 9 characters
  // after the decimal.
  return `_${Math.random()
    .toString(36)
    .substr(2, 9)}`;
};

const getStorageWithExpire = (cacheKey, expiry) => {
  const cached = localStorage.getItem(cacheKey);
  const whenCached = localStorage.getItem(`${cacheKey}:ts`);
  if (cached !== null && whenCached !== null) {
    // it was in sessionStorage! Yay!
    // Even though 'whenCached' is a string, this operation
    // works because the minus sign converts the
    // string to an integer and it will work.
    const age = (Date.now() - whenCached) / 1000;
    if (age < expiry) {
      return cached;
    }
    // We need to clean up this old key
    localStorage.removeItem(cacheKey);
    localStorage.removeItem(`${cacheKey}:ts`);
  }
  return null;
};

const getCachedSpectrum = async (id, options) => {
  let expiry = 10 * 60; // 10 min default
  if (typeof options === "number") {
    expiry = options;
  } else if (typeof options === "object") {
    // I hope you didn't set it to 0 seconds
    expiry = options.seconds || expiry;
  }

  const cacheKey = `spectrum_${id}`;
  const cached = getStorageWithExpire(cacheKey, expiry);
  if (cached !== null) return Promise.resolve(JSON.parse(cached));

  const { loading, error, data } = await client.query({
    query: GET_SPECTRUM,
    variables: { id }
  });
  if (!loading && !error) {
    const spec = data.spectrum;
    localStorage.setItem(cacheKey, JSON.stringify(spec));
    localStorage.setItem(`${cacheKey}:ts`, Date.now());
    return Promise.resolve(spec);
  }
};

async function fetchSpectra(id) {
  const data = await getCachedSpectrum(id);
  return {
    ...data,
    name: data.owner.name,
    inverted: false,
    visible: true,
    ecNormed: false,
    qyNormed: false
  };
}

const updateSpectra = async ({ toAdd, toRemove, dispatch }) => {
  let add = [];
  let remove = [];
  if (toAdd) {
    add = await Promise.all(toAdd.map(spec => fetchSpectra(spec.id)));
  }
  if (toRemove) {
    remove = toRemove.map(item => item.id);
  }
  dispatch({ type: "updateSeries", add, remove });
};

export { ID, getCachedSpectrum, getStorageWithExpire, updateSpectra };
