import gql from "graphql-tag"

const typeDefs = gql`
  extend type Query {
    activeSpectra: [Int]
    overlap: Spectrum
  }

  # extend type Mutation {
  #   setActiveSpectra(activeSpectra: [Int]!): [Int],
  #   updateActiveSpectra(add: [String], remove: [String]): [Int]
  # }

  extend type Spectrum {
    area: Float
  }

  type Selector {
    id: Int
    owner: String
    category: String
  }
`

export default typeDefs
