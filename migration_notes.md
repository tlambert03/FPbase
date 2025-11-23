# Migration notes

 - `common name` (VARCHAR 255)
 - `molecular formula` in Hill notation (VARCHAR)
 - `molecular weight` (DECIMAL 10,4)
 - `canonical SMILES` (TEXT) as the primary structure representation
 - `InChI` (TEXT) for standardized structure encoding
 - `InChIKey` (CHAR 27) provide uniqueness guarantees and cross-database compatibility
 - `CAS numbers` (VARCHAR 50) enable literature searches
 - `PubChem CIDs` (INTEGER) offer free cross-referencing to the world's largest chemical database.
