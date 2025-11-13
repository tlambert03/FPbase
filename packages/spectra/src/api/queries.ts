/**
 * GraphQL query strings
 * Replaces the old gql-tagged queries from Apollo
 */

// Fragment for fluorophore-specific fields
export const FLUOROPHORE_PARTS = `
  fragment FluorophoreParts on FluorophoreInterface {
    qy
    extCoeff
    twopPeakgm
    exMax
    emMax
  }
`

// Get a single spectrum by ID
export const GET_SPECTRUM = `
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
        ...FluorophoreParts
      }
    }
  }
  ${FLUOROPHORE_PARTS}
`

// Batch fetch multiple spectra
export function batchSpectraQuery(ids: string[]): string {
  const queries = ids
    .filter(Boolean)
    .map(
      (id) => `
    spectrum_${id}: spectrum(id: ${id}) {
      id
      data
      category
      color
      subtype
      owner {
        slug
        name
        id
        ...FluorophoreParts
      }
    }
  `
    )
    .join("\n")

  return `
    query BatchSpectra {
      ${queries}
    }
    ${FLUOROPHORE_PARTS}
  `
}

// Get optical configuration by ID
export const GET_OPTICAL_CONFIG = `
  query OpticalConfig($id: Int!) {
    opticalConfig(id: $id) {
      id
      name
      microscope {
        id
        name
      }
      filters {
        id
        path
        reflects
        spectrum {
          id
        }
      }
      light {
        id
        spectrum {
          id
        }
      }
      camera {
        id
        spectrum {
          id
        }
      }
      laser
      comments
    }
  }
`

// List all spectra (for search/autocomplete)
export const SPECTRA_LIST = `
  query SpectraList {
    spectra {
      id
      category
      subtype
      owner {
        id
        name
        slug
        url
      }
    }
  }
`

// TODO: autogenerate types from GraphQL schema in backend

interface SpectraSlug {
  id: string
  category: string
  subtype: string
  owner: {
    id: string
    name: string
    slug: string
    url: string | null
  }
}

export interface SpectraListResponse {
  spectra: SpectraSlug[]
}

// List optical configs (for search modal)
export const OPTICAL_CONFIG_LIST = `
  query OpticalConfigList {
    opticalConfigs {
      id
      name
      comments
      microscope {
        name
      }
    }
  }
`
