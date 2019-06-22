import gql from "graphql-tag"

const typeDefs = gql`
  extend type Query {
    activeSpectra: [Int]
  }
  extend type Mutation {
    setActiveSpectra(activeSpectra: [Int]!): [Int],
    updateActiveSpectra(add: [String], remove: [String]): [Int]
  }
`

export { typeDefs }
