import { library, dom } from '@fortawesome/fontawesome-svg-core';
import {
    faBook,
    faCamera,
    faCheck,
    faCog,
    faDownload,
    faExclamationCircle,
    faExternalLinkAlt,
    faEye,
    faFilter,
    faFlag,
    faHeart,
    faInfo,
    faInfoCircle,
    faLightbulb,
    faList,
    faMinusCircle,
    faPlus,
    faQuoteLeft,
    faSearch,
    faShareSquare,
    faTable,
    faTh,
    faThumbtack,
    faTimesCircle,
    faTrashAlt,
    faUpload,
    faWrench,
} from '@fortawesome/free-solid-svg-icons'; // ES Module "as" syntax
import { faHeart as farHeart, faEye as farEye } from '@fortawesome/free-regular-svg-icons';
import { faTwitter } from '@fortawesome/free-brands-svg-icons';

library.add(
    faBook,
    faCamera,
    faCheck,
    faCog,
    faDownload,
    faExclamationCircle,
    faExternalLinkAlt,
    faEye,
    faFilter,
    faFlag,
    faHeart,
    faInfo,
    faInfoCircle,
    faLightbulb,
    faList,
    faMinusCircle,
    faPlus,
    faQuoteLeft,
    faSearch,
    faShareSquare,
    faTable,
    faTh,
    faThumbtack,
    faTimesCircle,
    faTrashAlt,
    faUpload,
    faWrench,

    farHeart,
    farEye,

    faTwitter
);

// Replace any existing <i> tags with <svg> and set up a MutationObserver to
// continue doing this as the DOM changes.
dom.watch();
