import { gql } from "apollo-boost";

const GET_SPECTRUM = gql`
  query Spectrum($id: Int!) {
    spectrum(id: $id) {
      id
      data
      category
      color
      subtype
      owner {
        slug
        name
        id
      }
    }
  }
`;

const SPECTRA_LIST = gql`
  {
    spectra {
      id
      category
      subtype
      owner {
        name
        slug
        url
      }
    }
  }
`;

export { GET_SPECTRUM, SPECTRA_LIST };
