import * as React from 'react';
import { createRoot } from 'react-dom/client';
import { SearchAutocomplete } from '../components/SearchAutocomplete.jsx';

// Export function to initialize autocomplete on existing input
export function initAutocomplete() {
  // Find the search input (exists on both home page and navbar)
  const searchInput = document.querySelector('#algolia-search-input');
  if (!searchInput) {
    console.warn('Search input #algolia-search-input not found');
    return;
  }

  // Replace the input with autocomplete container
  const container = document.createElement('div');
  container.className = 'autocomplete-container';
  searchInput.parentElement.replaceChild(container, searchInput);

  // Initialize autocomplete
  const root = createRoot(container);
  root.render(React.createElement(SearchAutocomplete, { container }));
}
