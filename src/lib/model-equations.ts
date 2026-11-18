// These are the continuous-time equations whose unit-wise Euler updates are used by the simulators.
export type ModelEquation = { label: string; latex: string };

const dynamics: ModelEquation[] = [
  { label: 'Activity and state-specific adaptation', latex: String.raw`\tau_1\dot r_i^k=h_i^k-\theta_i^k-r_i^k,\qquad \tau_2\dot\theta_i^k=\sigma_i^k-\theta_i^k` },
  { label: 'Fast and slow inactive-state thresholds', latex: String.raw`T_{3A}\dot\theta_{i}^{0A}=\gamma_A(1-\sigma_i^0)-\theta_i^{0A},\qquad T_{3B}\dot\theta_{i}^{0B}=(1-\gamma_A)(1-\sigma_i^0)-\theta_i^{0B}` },
  { label: 'Potts probabilities', latex: String.raw`\sigma_i^k=\frac{e^{\beta r_i^k}}{Z_i},\quad \sigma_i^0=\frac{e^{\beta(U+\theta_i^{0A}+\theta_i^{0B})}}{Z_i},\quad Z_i=\sum_{l=1}^{S}e^{\beta r_i^l}+e^{\beta(U+\theta_i^{0A}+\theta_i^{0B})}` }
];

const localField = String.raw`h_i^k=\sum_{j,l}J_{ij}^{kl}\sigma_j^l+w\left(\sigma_i^k-\frac{1-\sigma_i^0}{S}\right)+I_i^k(t)`;
const regionalField = String.raw`h_{i,\rho}^k=\sum_{\nu\in\{F,P\}}\sum_{j,l}J_{ij,\rho\leftarrow\nu}^{kl}\sigma_{j,\nu}^l+w_\rho\left(\sigma_{i,\rho}^k-\frac{1-\sigma_{i,\rho}^0}{S_\rho}\right)+I_{i,\rho}^k(t)`;
const regionalCouplings = String.raw`J_{ij,\rho\leftarrow\rho}^{kl}=\frac{1+\lambda_\rho}{2C_\rho a_\rho(1-a_\rho/S_\rho)}\sum_\mu\eta_{i,\rho}^{\mu k}\eta_{j,\rho}^{\mu l}`;
const crossCouplings = String.raw`J_{ij,\rho\leftarrow\nu}^{kl}=\frac{1-\lambda_\rho}{2C_\rho\sqrt{a_\rho a_\nu(1-a_\rho/S_\rho)(1-a_\nu/S_\nu)}}\sum_{\mu,\kappa}A_{\mu\kappa}\eta_{i,\rho}^{\mu k}\eta_{j,\nu}^{\kappa l},\quad\nu\ne\rho`;

export const modelEquations = {
  homo: {
    rows: [
      { label: 'Recurrent field and cue', latex: localField },
      { label: 'Diluted Hebbian connections', latex: String.raw`J_{ij}^{kl}=\frac{1}{Ca(1-a/S)}\sum_\mu\eta_i^{\mu k}\eta_j^{\mu l},\qquad \eta_i^{\mu k}=\mathbf1[\xi_i^\mu=k]-\frac aS` },
      ...dynamics
    ],
    note: 'The connection sum uses C selected presynaptic units per target. k = 1,…,S denotes active states and 0 the inactive state. The cue I decays at the start of a run.'
  },
  hybrid: {
    rows: [
      { label: 'Field in each region', latex: regionalField },
      { label: 'Within-region Hebbian connections', latex: regionalCouplings },
      { label: 'Cross-region paired memories', latex: crossCouplings },
      ...dynamics
    ],
    note: 'F and P denote frontal and posterior regions. Here Aμκ = 1 only for a paired memory (μ = κ). η is the centred pattern indicator 1[ξ = k] − a/S. Each target samples C presynaptic units from each source region; λ belongs to the receiving region. Subscripts for region-specific dynamical parameters are omitted in the last three rows.'
  },
  topology: {
    rows: [
      { label: 'Field in each region', latex: regionalField },
      { label: 'Within-region Hebbian connections', latex: regionalCouplings },
      { label: 'Topology-weighted cross-region connections', latex: crossCouplings },
      { label: 'Association weights', latex: String.raw`A_{\mu\kappa}\geq0,\qquad\sum_\kappa A_{\mu\kappa}=1` },
      ...dynamics
    ],
    note: 'A is the selected reciprocal memory association matrix; the original topology has one-to-one pairs. η is the centred pattern indicator 1[ξ = k] − a/S. Only sampled presynaptic connections contribute, and λ belongs to the receiving region. Subscripts for region-specific dynamical parameters are omitted in the last three rows.'
  }
} as const;
