#!/usr/bin/env node
/**
 * Extract Font Awesome icons used in FPbase codebase
 * Generates a TypeScript file with SVG path data for all icons
 */

import { writeFileSync } from 'fs';
import { join } from 'path';

// Icon mapping: kebab-case class names → camelCase FA exports
const iconMap = {
  // Solid icons
  'arrow-left': { pkg: 'solid', name: 'faArrowLeft' },
  'arrow-right': { pkg: 'solid', name: 'faArrowRight' },
  'arrows': { pkg: 'solid', name: 'faArrows' },
  'bars': { pkg: 'solid', name: 'faBars' },
  'book': { pkg: 'solid', name: 'faBook' },
  'camera': { pkg: 'solid', name: 'faCamera' },
  'chart-area': { pkg: 'solid', name: 'faChartArea' },
  'check': { pkg: 'solid', name: 'faCheck' },
  'check-circle': { pkg: 'solid', name: 'faCheckCircle' },
  'clock': { pkg: 'solid', name: 'faClock' },
  'code': { pkg: 'solid', name: 'faCode' },
  'cog': { pkg: 'solid', name: 'faCog' },
  'dot-circle-o': { pkg: 'regular', name: 'faDotCircle' }, // v4 → v5
  'download': { pkg: 'solid', name: 'faDownload' },
  'edit': { pkg: 'solid', name: 'faEdit' },
  'envelope': { pkg: 'solid', name: 'faEnvelope' },
  'exchange': { pkg: 'solid', name: 'faExchange' },
  'exchange-alt': { pkg: 'solid', name: 'faExchangeAlt' },
  'exclamation-circle': { pkg: 'solid', name: 'faExclamationCircle' },
  'exclamation-triangle': { pkg: 'solid', name: 'faExclamationTriangle' },
  'external-link-alt': { pkg: 'solid', name: 'faExternalLinkAlt' },
  'eye': { pkg: 'solid', name: 'faEye' },
  'fast-backward': { pkg: 'solid', name: 'faFastBackward' },
  'fast-forward': { pkg: 'solid', name: 'faFastForward' },
  'filter': { pkg: 'solid', name: 'faFilter' },
  'flag': { pkg: 'solid', name: 'faFlag' },
  'flash': { pkg: 'solid', name: 'faBolt' }, // v4 → v5 renamed
  'flip-horizontal': { pkg: 'solid', name: 'faArrowsLeftRight' }, // closest match
  'floppy-o': { pkg: 'solid', name: 'faFloppyDisk' }, // v4 → v5
  'heart': { pkg: 'solid', name: 'faHeart' },
  'heart-o': { pkg: 'regular', name: 'faHeart', exportName: 'farHeart' }, // outline heart for favorites
  'home': { pkg: 'solid', name: 'faHome' },
  'info': { pkg: 'solid', name: 'faInfo' },
  'info-circle': { pkg: 'solid', name: 'faInfoCircle' },
  'key': { pkg: 'solid', name: 'faKey' },
  'keyboard': { pkg: 'solid', name: 'faKeyboard' },
  'lightbulb': { pkg: 'solid', name: 'faLightbulb' },
  'link': { pkg: 'solid', name: 'faLink' },
  'list': { pkg: 'solid', name: 'faList' },
  'minus': { pkg: 'solid', name: 'faMinus' },
  'minus-circle': { pkg: 'solid', name: 'faMinusCircle' },
  'pause': { pkg: 'solid', name: 'faPause' },
  'play': { pkg: 'solid', name: 'faPlay' },
  'plus': { pkg: 'solid', name: 'faPlus' },
  'plus-circle': { pkg: 'solid', name: 'faPlusCircle' },
  'power-off': { pkg: 'solid', name: 'faPowerOff' },
  'question-circle': { pkg: 'solid', name: 'faQuestionCircle' },
  'quote-left': { pkg: 'solid', name: 'faQuoteLeft' },
  'search': { pkg: 'solid', name: 'faSearch' },
  'share': { pkg: 'solid', name: 'faShare' },
  'share-square': { pkg: 'solid', name: 'faShareSquare' },
  'spinner': { pkg: 'solid', name: 'faSpinner' },
  'square': { pkg: 'solid', name: 'faSquare' },
  'square-o': { pkg: 'regular', name: 'faSquare' }, // v4 → v5 regular
  'star': { pkg: 'solid', name: 'faStar' },
  'step-backward': { pkg: 'solid', name: 'faStepBackward' },
  'step-forward': { pkg: 'solid', name: 'faStepForward' },
  'sun': { pkg: 'solid', name: 'faSun' },
  'table': { pkg: 'solid', name: 'faTable' },
  'th': { pkg: 'solid', name: 'faTh' },
  'thumbtack': { pkg: 'solid', name: 'faThumbtack' },
  'times': { pkg: 'solid', name: 'faTimes' },
  'times-circle': { pkg: 'solid', name: 'faTimesCircle' },
  'trash': { pkg: 'solid', name: 'faTrash' },
  'trash-alt': { pkg: 'solid', name: 'faTrashAlt' },
  'undo': { pkg: 'solid', name: 'faUndo' },
  'upload': { pkg: 'solid', name: 'faUpload' },
  'user': { pkg: 'solid', name: 'faUser' },
  'wrench': { pkg: 'solid', name: 'faWrench' },

  // Brand icons
  'github': { pkg: 'brands', name: 'faGithub' },
  'google': { pkg: 'brands', name: 'faGoogle' },
  'orcid': { pkg: 'brands', name: 'faOrcid' },
  'x-twitter': { pkg: 'brands', name: 'faXTwitter' },
};

async function extractIcons() {
  const icons = [];
  const errors = [];

  for (const [kebabName, config] of Object.entries(iconMap)) {
    const { pkg, name, exportName } = config;
    try {
      const pkgName = `@fortawesome/free-${pkg}-svg-icons`;
      const module = await import(pkgName);
      const icon = module[name];

      if (!icon) {
        errors.push(`Icon ${name} not found in ${pkgName}`);
        continue;
      }

      const [width, height, , , svgPathData] = icon.icon;

      // Handle both single path and array of paths
      const path = Array.isArray(svgPathData) ? svgPathData.join(' ') : svgPathData;

      // Generate consistent prefix based on package: fas/far/fab
      const prefix = pkg === 'solid' ? 'fas' : pkg === 'regular' ? 'far' : 'fab';
      const consistentName = exportName || name.replace(/^fa/, prefix);

      icons.push({
        kebabName,
        exportName: consistentName,
        width,
        height,
        path,
      });
    } catch (err) {
      errors.push(`Error importing ${name}: ${err.message}`);
    }
  }

  if (errors.length > 0) {
    console.error('Errors encountered:');
    errors.forEach(err => console.error(`  - ${err}`));
  }

  return icons;
}

function generateTypeScript(icons) {
  const header = `/**
 * Static FontAwesome icon data
 * Extracted from @fortawesome packages to eliminate runtime dependency
 *
 * Icon data format: { width, height, path }
 * - width/height: SVG viewBox dimensions
 * - path: SVG path data string
 *
 * Generated by scripts/extract-fa-icons.mjs
 */

export interface IconData {
  width: number
  height: number
  path: string
}
`;

  const iconExports = icons
    .sort((a, b) => a.exportName.localeCompare(b.exportName))
    .map(({ exportName, width, height, path }) => {
      // Escape quotes in path data
      const escapedPath = path.replace(/\\/g, '\\\\').replace(/"/g, '\\"');

      return `export const ${exportName}: IconData = {
  width: ${width},
  height: ${height},
  path: "${escapedPath}",
}`;
    })
    .join('\n\n');

  return header + '\n' + iconExports + '\n';
}

function generatePython(icons) {
  const header = `"""
Static FontAwesome icon data
Extracted from @fortawesome packages to eliminate runtime dependency

Icon data format: dict with 'width', 'height', 'path' keys
- width/height: SVG viewBox dimensions
- path: SVG path data string

Generated by scripts/extract-fa-icons.mjs
"""

from __future__ import annotations

from typing import TypedDict


class IconData(TypedDict):
    """SVG icon data."""

    width: int
    height: int
    path: str


`;

  // Create dict entries
  const iconDict = icons
    .sort((a, b) => a.exportName.localeCompare(b.exportName))
    .map(({ exportName, kebabName, width, height, path }) => {
      // Escape quotes and backslashes for Python
      const escapedPath = path.replace(/\\/g, '\\\\').replace(/"/g, '\\"');

      return `    "${exportName}": {
        "width": ${width},
        "height": ${height},
        "path": "${escapedPath}",
    }`;
    })
    .join(',\n');

  // Also create a kebab-case alias mapping for template tag usage
  const aliasMapping = icons
    .sort((a, b) => a.kebabName.localeCompare(b.kebabName))
    .map(({ exportName, kebabName }) => `    "${kebabName}": "${exportName}"`)
    .join(',\n');

  return `${header}
ICONS: dict[str, IconData] = {
${iconDict}
}

# Kebab-case aliases for template tag convenience
# Maps 'trash' → 'fasTrash', 'square-o' → 'farSquare', etc.
ICON_ALIASES: dict[str, str] = {
${aliasMapping}
}


def get_icon(name: str) -> IconData | None:
    """
    Get icon data by name.

    Accepts both camelCase exports (fasTrash) and kebab-case aliases (trash).
    """
    # Try direct lookup first
    if name in ICONS:
        return ICONS[name]

    # Try alias lookup
    if name in ICON_ALIASES:
        return ICONS[ICON_ALIASES[name]]

    return None
`;
}

async function main() {
  console.log('Extracting Font Awesome icons...');
  const icons = await extractIcons();
  console.log(`✓ Extracted ${icons.length} icons`);

  // Generate TypeScript for frontend
  const tsContent = generateTypeScript(icons);
  const tsPath = join(process.cwd(), 'frontend/src/icons/fa-icons.ts');
  writeFileSync(tsPath, tsContent, 'utf8');
  console.log(`✓ Written TypeScript to ${tsPath}`);

  // Generate Python for Django templates
  const pyContent = generatePython(icons);
  const pyPath = join(process.cwd(), 'backend/fpbase/icons.py');
  writeFileSync(pyPath, pyContent, 'utf8');
  console.log(`✓ Written Python to ${pyPath}`);
}

main().catch(err => {
  console.error('Fatal error:', err);
  process.exit(1);
});
